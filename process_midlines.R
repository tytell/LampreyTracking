require(tidyverse)

process_filepath <- function(filepath) {
  df <- str_match(filepath, 'Animal(\\d)/(Control|\\d[dh]pi)/(\\d{6})')
  colnames(df) <- c("all", "Animal", "Treatment", "Date")
  df %>%
    as_tibble() %>%
    mutate(Date = lubridate::mdy(Date),
           Treatment = factor(Treatment, levels = treatments)) %>%
    select(-all)
}

load_katz_data <- function(filename, widthdata, fps, bodypartorder, showfile=FALSE) {
  if (showfile) {
    print(basename(filename))
  }
  df <- read_csv(filename, show_col_types = FALSE) %>%
    bind_cols(process_filepath(filename)) %>%
    left_join(scales, by = c('Animal', 'Treatment', 'Date')) %>%
    mutate(t = frame / fps,
           xmm = x * Scale,
           ymm = y * Scale,
           bodyparts = factor(bodyparts, levels = bodypartorder),
           bodypartord = as.numeric(bodyparts),
           fullpathname = filename,
           filename = basename(filename),
           trial = str_extract(filename, '\\d{6}_A\\d+V\\d+_\\d+'))
}

get_arc_length <- function(df) {
  #' Calculates arc length along the body
  
  df %>%
    group_by(frame, .add=TRUE) %>%
    mutate(dx = lead(xmm) - xmm,
           dy = lead(ymm) - ymm,
           s = cumsum(sqrt(dx^2 + dy^2)),
           s = replace_na(lag(s), 0)) %>%
    ungroup()
}

get_median_length <- function(df) {
  #' Pulls out the last value of arc length in each frame and uses that to
  #' estimate the median body length across all frames

  df %>%
    group_by(frame, .add=TRUE) %>%
    summarize(len = last(s)) %>%
    summarize(len = median(len)) %>%
    pull(len)
}

interpolate_width <- function(df, widthdata) {
  #' Interpolates the width of the animal
  #' 
  #' @param df Data frame containing s, xmm, ymm, the arc length along the body, 
  #' and coordinates along the body in mm. Points do not need to be evenly spaced.
  #' @param widthdata Vector containing the width of the body at evenly spaced
  #' positions along the body, normalized to body length
  
  s0 <- pracma::linspace(0, 1, length(widthdata))
  widthfun <- approxfun(s0, widthdata, yleft=0, yright=0)
  
  len <- get_median_length(df)
  
  df %>%
    mutate(width = widthfun(s/len) * len)
}

get_center_of_mass <- function(df) {
  #' Computes the center of mass, assuming constant density
  #' 
  #' @param df Kinematics data frame
  
  height <- max(df$width)
  
  com <-
    df %>%
    group_by(frame, .add=TRUE) %>%
    summarize(m = pracma::trapz(s, pi*height*width/4),
              comx = pracma::trapz(s, xmm * pi*height*width/4) / m,
              comy = pracma::trapz(s, ymm * pi*height*width/4) / m) %>%
    ungroup()
  
  df %>%
    left_join(com, by = "frame")
}

smooth_com_spline <- function(t,com, spar) {
  #' Smooths the center of mass location using a smoothing spline
  #' 
  #' @param t Time
  #' @param com Location of COM (in x or y)
  #' @param spar Smoothing parameter. Seems to work well with a value of 0.8
  
  coms <- numeric(length(com))
  good <- !is.na(com)
  
  sp <- smooth.spline(t[good], com[good], spar = spar)
  
  coms[!good] <- NA
  coms[good] = predict(sp)$y
  
  coms
}

get_swim_vel_dir <- function(df, s = 0.8) {
  #' Computes the swimming velocity and a smoothed swimming direction vector.
  #' The swimming direction vector is used later on for calculating amplitudes
  #' and excursions.
  #' 
  #' @param df Data frame
  #' @param s Smoothing parameter. Seems to work will around 0.8
  #' @param fps Frames per second
  
  swim <-
    df %>%
    group_by(frame, .add=TRUE) %>%
    # COM has one location per frame, but we have multiple body positions, so we
    # just take the first value
    summarize(across(c(comx, comy, t), first)) %>%
    mutate(swimvelx = 0.5 * ((lead(comx) - comx) / (lead(t) - t) + (comx - lag(comx)) / (t - lag(t))),   # central difference derivative
           swimvely = 0.5 * ((lead(comy) - comy) / (lead(t) - t) + (comy - lag(comy)) / (t - lag(t))),
           swimvel = sqrt(swimvelx^2 + swimvely^2),  # magnitude is speed
           swimvelxs = smooth_com_spline(t, swimvelx, s),    # smooth the COM x location
           swimvelys = smooth_com_spline(t, swimvely, s),
           swimvels = sqrt(swimvelxs^2 + swimvelys^2),
           swimdirx = swimvelxs / swimvels,   # make a normal vector corresponding to the swimming direction
           swimdiry = swimvelys / swimvels) %>%
    select(-comx, -comy, -t)
  
  # merge back up with the main data frame
  df %>%
    left_join(swim, by = c("frame"))
}

get_excursions <- function(df) {
  #' Gets lateral excursion of the body, relative to the swimming direction.
  #' 
  #' Projects the position of each body segment, relative to the center of mass, on
  #' to the vector perpendicular to the swimming direction
  df %>%
    mutate(exc = -(xmm - comx) * swimdiry + 
             (ymm - comy) * swimdirx)
}

get_cycles <- function(df) {
  #' Uses the time series for lateral excursion to find the cycle periods
  #' 
  #' For each body segment, it looks for where the excursion crosses the swimming
  #' direction (a zero crossing). Uses linear interpolation to estimate the actual
  #' time of the zero crossing.
  #' 
  #' Then uses those times to estimate the period and frequency of the movement
  #' at each body segment.
  #' 
  #' Finally, estimates a cycle phase for for the tail, which we use as an
  #' overall phase estimate.
  df <-
    df %>%
    group_by(bodyparts, .add=TRUE) %>%
    mutate(excn = lead(exc),
           zerocross = case_when(sign(excn) > sign(exc)   ~   1,
                                 sign(excn) == sign(exc)  ~   0,
                                 sign(excn) < sign(exc)   ~   -1,
                                 TRUE   ~   0)) %>%
    ungroup(bodyparts)
  
  zerocross <-
    df %>%
    group_by(bodyparts, .add=TRUE) %>%
    filter(zerocross != 0) %>%
    mutate(t0 = t + (1/fps)/(excn - exc) * (0-exc),
           per = lead(t0)-lag(t0),
           freqbody = 1/per) %>%
    select(bodyparts, frame, t0, per,freqbody)

  df <-
    df %>%
    left_join(zerocross, by = c("frame", "bodyparts"))
  
  t0tail <-
    df %>%
    ungroup() %>%
    filter(bodyparts == "tailtip" &
             zerocross != 0) %>%
    select(t0, freqbody) %>%
    mutate(phtail = seq(from = 0, length.out = n())/2)

  df %>%
    ungroup() %>%
    mutate(freq = approx(t0tail$t0, t0tail$freqbody, t)$y,
           phase = approx(t0tail$t0, t0tail$phtail, t)$y,
           cycle = floor(phase*2) / 2)
}

get_amplitudes <- function(df) {
  #' Computes body undulation amplitude.
  #' 
  #' Saves the maximum absolute amplitude in each half cycle and the
  #' frame at which it occurs. Then smooths the amplitude using the 
  #' two amplitude estimates on either side.
  amps <-
    df %>%
    group_by(bodyparts, cycle) %>%
    summarize(ampframe = frame[which.max(abs(exc))],
              amp = max(abs(exc))) %>%
    ungroup() %>%
    mutate(amp = (lag(amp) + 2*amp + lead(amp)) / 4)

  left_join(df, amps, by = c("bodyparts", "cycle", "frame" = "ampframe"))
}

get_next_level <- function(fct) {
  #' Gets the next level in an ordered factor.
  #' 
  #' Returns NA if you try to get the next level after the last one
  lev <- levels(fct)
  n = as.numeric(fct) + 1
  
  good = (n >= 1) & (n <= length(lev))
  n[!good] = 1
  
  nextfct <- factor(lev[n], levels = lev)
  nextfct[!good] = NA
  nextfct
}

get_wavespeed <- function(df) {
  #' Estimates the wavespeed from the swimming kinematics.
  #' 
  #' Looks at each segment and chooses the time of the current zero crossing
  #' and the next zero crossing. Then looks at the next segment along the body and
  #' finds the zero crossing in that segment that has the same sign as the current
  #' zero crossing. The wave speed is the arc length between the two segments
  #' divided by the difference in time of the zero crossings.
 
  bodyparts <- levels(df$bodyparts)
  bodyparts <- factor(bodyparts, levels = bodyparts)

  # run through the segments
  t0seg = list()
  for (i in seq(length.out = length(bodyparts)-1)) {
    seg1 <- bodyparts[i]
    
    seg2 = get_next_level(seg1)
    
    # for this segment, save the time of the next zero crossing
    # as a new column
    t0seg1 <- df %>%
      ungroup() %>%
      filter(bodyparts == seg1 & !is.na(t0)) %>%
      arrange(t0) %>%
      mutate(t0next = lead(t0)) %>%
      select(bodyparts, s, t0, t0next, zerocross)
    
    # for the next segment, save the arc length of the segment and
    # the times of the zero crossings
    t0seg2 = df %>%
      ungroup() %>%
      filter(bodyparts == seg2 & !is.na(t0)) %>%
      select(bodyparts, s, t0, zerocross)
    
    # then run through all of the zero crossings in the current segment
    t0seg1$t0nextseg <- NA
    t0seg1$snextseg <- NA
    for (r in seq(length.out = nrow(t0seg1))) {
      # in the next segment, look for a zero crossing that's after the one
      # in the current segment, but before the next one in the current segment, and
      # has the same sign zero crossing.
      nextt0 <- 
        t0seg2 %>%
        filter((t0 >= t0seg1$t0[r]) & (t0 < t0seg1$t0next[r]) & (zerocross == t0seg1$zerocross[r])) %>%
        transmute(t0nextseg = t0,
                  snextseg = s) %>%
        filter(!is.na(t0nextseg))
      

      # save out the times of the next segment zero crossing and the arc length of the next segment
      t0seg1$t0nextseg[r] <- pull(nextt0, t0nextseg) %>% first()
      t0seg1$snextseg[r] <- pull(nextt0, snextseg) %>% first()
    }
    
    # also save the name of the next segment
    t0seg1$nextseg = seg2
    t0seg[[i]] <- t0seg1 %>%
      mutate(wavespeed = (snextseg - s)/(t0next - t0))
  }
  
  # stack all of the zero crossing data
  t0seg <- bind_rows(t0seg)
  
  # and merge it in with the main data set
  df %>%
    left_join(t0seg, by = c("bodyparts", "t0", "zerocross", "s"))
}

get_wavelength <- function(df) {
  #' Get the wavelength along the body.
  #' 
  #' For each frame, look for nodes (points along the body where the excursion
  #' goes from one sign to the other). We have relatively few segments, so it
  #' uses a linear interpolation to find the actual location of the node.
  #' 
  #' The wavelength is twice the difference in arc length between two nodes of opposite
  #' sign along the body. (or just the difference in arc length between two nodes of the
  #' same sign, but it's rare that we would have that)
  df <-
    df %>%
    group_by(frame) %>%
    arrange(bodyparts, .by_group = TRUE) %>%
    mutate(excnextseg = lead(exc),
           snextseg = lead(s),
           node = case_when(sign(lead(exc)) > sign(exc)   ~   1,
                            sign(lead(exc)) < sign(exc)   ~   -1,
                            TRUE   ~   0))
  
  nodes <-
    df %>%
    group_by(frame) %>%
    filter(node != 0) %>%
    mutate(s0 = s + (snextseg - s) / (excnextseg - exc) * (0 - exc),
           wavelength = lead(s0) - s0,
           wavelength = case_when(sign(lead(node)) != sign(node)   ~   wavelength*2,
                                  sign(lead(node)) == sign(node)   ~   wavelength)) %>%
    select(frame, bodyparts, s0, wavelength)
  
  df %>%
    left_join(nodes, by = c("frame", "bodyparts"))
}

get_all_kinematics <- function(df,
                               widthdata) {
  df %>%
    get_arc_length() %>%
    interpolate_width(widthdata) %>%
    get_center_of_mass() %>%
    get_swim_vel_dir() %>%
    get_excursions() %>%
    get_cycles() %>%
    get_amplitudes() %>%
    get_wavespeed() %>%
    get_wavelength()
}

get_point_order <- function(df) {
  df %>%
    group_by(trial, frame) %>%
    group_modify(~ get_point_order_in_frame(.x))
}

get_point_order_in_frame <- function(df) {
  xy <- select(df, xmm, ymm)
  D <- dist(xy)
  D <- as.matrix(D)
  
  n <- dim(D)[1]
  ind <- rep(NA_integer_, n)
  for (i in seq_len(n-1)) {
    d1 <- D[(i+1):n,i]
    ind[i] = which.min(d1) + i
  }
  
  df$nextpoint <- ind
  
  df
}
