require(tidyverse)
require(RcppRoll)

process_filepath <- function(filepath) {
  df <- str_match(filepath, 'Exp_(\\d{6})/Animal\\s*(\\d)/(Control|\\d[dh]pi)/(\\d{6})')
  colnames(df) <- c("all", "Date", "Animal", "Treatment", "VideoDate")
  df %>%
    as_tibble() %>%
    mutate(Date = lubridate::mdy(Date),
           Treatment = factor(Treatment, levels = treatments),
           VideoDate = lubridate::mdy(VideoDate)) %>%
    select(-all)
}

load_katz_data <- function(filename, scales, widthdata, fps, bodypartorder, showfile=FALSE) {
  if (showfile) {
    print(basename(filename))
  }
  df <- read_csv(filename, show_col_types = FALSE) %>%
    bind_cols(process_filepath(filename)) 
  
  df %>%
    left_join(scales, by = c('Date', 'Animal', 'Treatment')) %>%
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

smooth_point_spline <- function(t,com, spar) {
  #' Smooths the center of mass location using a smoothing spline
  #' 
  #' @param t Time
  #' @param com Location of COM (in x or y)
  #' @param spar Smoothing parameter. Seems to work well with a value of 0.8
  
  coms <- numeric(length(com))
  good <- !is.na(com)
  
  if (any(good)) {
    sp <- smooth.spline(t[good], com[good], spar = spar)
    
    coms[good] = predict(sp)$y
  }  

  coms[!good] <- NA
  coms
}

get_swim_vel_dir <- function(df, s = 0.8) {
  #' Computes the swimming velocity and a smoothed swimming direction vector.
  #' The swimming direction vector is used later on for calculating amplitudes
  #' and excursions.
  #' 
  #' @param df Data frame
  #' @param s Smoothing parameter. Seems to work well around 0.8
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
           swimvelxs = smooth_point_spline(t, swimvelx, s),    # smooth the COM x location
           swimvelys = smooth_point_spline(t, swimvely, s),
           swimvels = sqrt(swimvelxs^2 + swimvelys^2),
           swimdirx = swimvelxs / swimvels,   # make a normal vector corresponding to the swimming direction
           swimdiry = swimvelys / swimvels) %>%
    select(-comx, -comy, -t)
  
  # merge back up with the main data frame
  df %>%
    left_join(swim, by = c("frame"))
}

get_body_axis_in_frame <- function(df) {
  #' Gets a central axis for one set of points in a frame.
  #' 
  #' Uses `svd` to estimate the central x and y axis.
  #' Requires `comx` and `comy` to already have been estimated in
  #' the data frame.
  
  XY <- cbind(df$xmm - df$comx, df$ymm - df$comy)
  svdout <- svd(XY)

  if (is.list(svdout)) {
    df <-
      df %>%
      mutate(bodyaxisx = svdout$v[1,1],
             bodyaxisy = svdout$v[2,1])
  } else {
    df <-
      df %>%
      mutate(bodyaxisx = NA,
             bodyaxisy = NA)
  }
  df
}

get_body_axis <- function(df, s = 0.8) {
  #' Gets the central axis of the body
  #' 
  #' Uses `get_body_axis_in_frame` to estimate the central body axis.
  #' Requires `comx` and `comy` to have already been estimated.
  #' 
  #' @param s Smoothing parameter for smoothing body axis fluctuations
  
  ax <-
    df %>%
    group_by(frame) %>%
    group_modify(~ get_body_axis_in_frame(.x))
  
  ax <-
    ax %>%
    group_by(frame) %>%
    summarize(across(c(t, bodyaxisx, bodyaxisy), first)) %>%
    mutate(signflip = if_else(bodyaxisx * lag(bodyaxisx) + bodyaxisy * lag(bodyaxisy) < -0.8, -1, 1),
           signflip = replace_na(signflip, 1))
  
  ax <-
    ax %>%
    ungroup() %>%
    mutate(signflip = cumprod(signflip),
           bodyaxisx = bodyaxisx * signflip,
           bodyaxisy = bodyaxisy * signflip)
    
  ax <-
    ax %>%
    mutate(bodyaxisxfr = bodyaxisx,
           bodyaxisyfr = bodyaxisy,
           bodyaxisx = smooth_point_spline(t, bodyaxisxfr, s),
           bodyaxisy = smooth_point_spline(t, bodyaxisyfr, s),
           bodyaxismag = sqrt(bodyaxisx^2 + bodyaxisy^2),
           bodyaxisx = bodyaxisx / bodyaxismag,
           bodyaxisy = bodyaxisy / bodyaxismag) %>%
    select(-t, -bodyaxismag)
    
  df %>%
    left_join(ax, by = c("frame"))
}

get_excursions <- function(df) {
  #' Gets lateral excursion of the body, relative to the swimming direction.
  #' 
  #' Projects the position of each body segment, relative to the center of mass, on
  #' to the vector perpendicular to the swimming direction
  df %>%
    mutate(exc = -(xmm - comx) * bodyaxisy + 
             (ymm - comy) * bodyaxisx)
}

unwrap <- function(y, modulus = 2*pi, jump = 0.5, na.rm = FALSE) {
  #' Unwraps angle data to remove discontinuities.
  #' 
  #' With angle data mod 2pi, the angles can jump discontinuously from
  #' pi to -pi or -pi to pi. This removes those discontinuities by adding
  #' or subtracting 2pi when there is a jump.
  #' 
  #' @param y Time series angle data.
  #' @param modulus Modulus of the data. For angle data in radians, use
  #'    2pi. Could also be 360 for angle data in degrees or 1 for data
  #'    mod 1.
  #' @param jump What fraction of the modulus constitutes a jump. Default
  #'    is 0.5, but you could set it up to 0.9 or more to be stricter.
  #' @param na.rm Remove NA values before calculating jumps (or not)
  
  if (na.rm) {
    good = !is.na(y)
  } else {
    good = rep_len(TRUE, length(y))
  }
  
  y1 <- y1[good]
  
  dy <- y1 - lag(y1)
  dy[1] <- 0
  
  step <- rep_len(0, length(y))
  
  # look for the jumps
  step[good] = case_when(dy < -jump*modulus  ~  modulus,
                   dy > jmup*modulus  ~  -modulus,
                   TRUE  ~  0)
  # steps are cumulative, so make sure to add them up
  step = cumsum(step)
  
  # then add the steps back on to the original data
  y + step
}

deriv <- function(x, y, ord = 1, method = 'direct') {
  #' Estimate first or second derivatives for dy/dx.
  #' 
  #' Uses central differencing where possible.
  #' 
  #' @param x x variable. Does not need to be evenly spaced.
  #' @param y y variable.
  #' @param ord Order of the derivative (1 or 2).
  #' @param method Method for taking second derivatives. Either
  #'   * 'direct' (default) Uses a direct formula, based on a central difference of
  #'     forward and backward differences, from [https://mathformeremortals.wordpress.com/2013/01/12/a-numerical-second-derivative-from-three-points/]
  #'   * 'repeat' Repeat two first derivatives.
  
  if (ord == 1) {
    # standard central difference formula for first derivative
    D <- (lead(y) - lag(y)) / (lead(x) - lag(x))
    
    # forward difference for the first point
    D[1] <- (y[2] - y[1]) / (x[2] - x[1])
    
    # backward difference for the second
    n <- length(x)
    D[n] <- (y[n] - y[n-1]) / (x[n] - x[n-1])
  } else if (ord == 2) {
    if (method == 'direct') {
      # direct formula for the second derivative, given uneven spacing in x
      # see https://mathformeremortals.wordpress.com/2013/01/12/a-numerical-second-derivative-from-three-points/
      D <- 2*lag(y) / ((x - lag(x))*(lead(x) - lag(x))) -
        2*y / ((lead(x) - x)*(x - lag(x))) +
        2*lead(y) / ((lead(x) - x)*(lead(x) - lag(x)))
    } else if (method == 'repeat') {
      # second derivative by repeating first derivatives
      dydx <- deriv(x, y, ord = 1)
      D <- deriv(x, dydx, ord = 1)
    }
  }
  D
}

get_curvature_in_frame <- function(s, x, y) {
  #' Estimates curvature for a single curve
  #' 
  #' Assumes that points are in order
  
  dx <- deriv(s, x)
  dy <- deriv(s, y)
  ddx <- deriv(s, x, 2, method='direct')
  ddy <- deriv(s, y, 2, method='direct')
  
  (dx*ddy - dy*ddx) / ((dx^2 + dy^2)^1.5)
}

get_curvature <- function(df) {
  df %>%
    arrange(frame, bodyparts) %>%
    group_by(frame) %>%
    mutate(curve = get_curvature_in_frame(s, xmm, ymm) * len)
}

get_cycles <- function(df, track = 'excursion',
                       include.zeros = FALSE,
                       smooth.excursion = 0.2,
                       min.peak.gap = 0.01,
                       min.peak.size = 0.01) {
  #' Uses the time series for curvature or lateral excursion to find the cycle periods
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
  
  # first check whether we're tracking curvature or excursion
  if (track == 'curvature') {
    df <-
      df %>%
      mutate(trackvar = curve)
  }
  else if (track == 'excursion') {
    df <-
      df %>%
      mutate(trackvar = exc)
  }
  else {
    stop(sprintf('Unknown tracking option %s', track))
  }
  
  # first smooth the excursion, then look for peaks and zerocrossings
  df <-
    df %>%
    group_by(bodyparts) %>%
    mutate(tracks = smooth_point_spline(t, trackvar, smooth.excursion),
           peak = case_when((tracks > lag(tracks)) & (tracks > lead(tracks))   ~  1,
                            (tracks < lag(tracks)) & (tracks < lead(tracks))   ~  -1,
                            TRUE   ~  0),
           zerocross = case_when(sign(lead(exc)) > sign(exc)   ~   1,
                                 sign(lead(exc)) < sign(exc)   ~   -1,
                                 TRUE   ~   0)) %>%
    ungroup()

  # now interpolate the exact time of the peaks and the zero crossings
  # assume the peak is a parabola
  # formula from https://stackoverflow.com/questions/717762/how-to-calculate-the-vertex-of-a-parabola-given-three-points   
  # assuming the zero crossing is a line
  df <-
    df %>%
    ungroup() %>%
    arrange(bodyparts, t) %>%
    group_by(bodyparts) %>%
    mutate(denom = if_else(peak != 0, 
                           (lag(t) - t) * (lag(t) - lead(t)) * (t - lead(t)), 
                           NA_real_),
           A = if_else(peak != 0, 
                       (lead(t) * (tracks - lag(tracks)) + t * (lag(tracks) - lead(tracks)) + lag(t) * (lead(tracks) - tracks)) / denom, 
                       NA_real_),
           B = if_else(peak != 0,
                       (lead(t)^2 * (lag(tracks) - tracks) + t^2 * (lead(tracks) - lag(tracks)) + lag(t)^2 * (tracks - lead(tracks))) / denom,
                       NA_real_),
           tp = -B / (2*A),
           t0 = if_else(zerocross != 0,
                        t + (lead(t) - t) / (lead(tracks) - tracks) * (0-tracks),
                        NA_real_)) %>%
    select(-denom, -A, -B)
  
  # get rid of peaks or zeros that are too close to one another
  cycleorder <-
    df %>%
    filter((peak != 0) | (zerocross != 0)) %>%
    mutate(cycle_type = case_when(peak != 0  ~  'peak',
                                  zerocross != 0  ~  'zero')) %>%
    group_by(bodyparts, cycle_type) %>%
    mutate(peakgapfwd = lead(tp) - tp,
           peakgaprev = tp - lag(tp),
           peakgap = case_when(is.na(peakgapfwd)  ~  peakgaprev,
                               is.na(peakgaprev)  ~  peakgapfwd,
                               peakgapfwd < peakgaprev  ~  peakgapfwd,
                               TRUE  ~  peakgaprev),
           tp = if_else(peakgap < min.peak.gap, NA_real_, tp))
  
  cycleorder <-
    cycleorder %>%
    mutate(peak = if_else(is.na(tp), 0, peak))
  
  # also get rid of peaks that are too small
  cycleorder <-
    cycleorder %>%
    filter((peak != 0) | (zerocross != 0)) %>%
    mutate(cycle_type = case_when(peak != 0  ~  'peak',
                                  zerocross != 0  ~  'zero')) %>%
    arrange(bodyparts, cycle_type, frame) %>%
    group_by(bodyparts, cycle_type) %>%
    mutate(peaksizerev = abs(exc - lag(exc)),
           peaksizefwd = abs(exc - lead(exc)),
           # take the smaller of the forward and backward differences
           peaksize = case_when(is.na(peaksizerev)  ~  peaksizefwd,
                                is.na(peaksizefwd)  ~  peaksizerev,
                                peaksizerev > peaksizefwd  ~  peaksizefwd,
                                TRUE  ~  peaksizerev),
           peaksize = if_else(cycle_type == 'zero', NA_real_, peaksize) / len)
  
  cycleorder <-
    cycleorder %>%
    group_by(bodyparts) %>%
    mutate(peak = if_else((cycle_type == 'peak') & (peaksize < min.peak.size), 0, peak))    

  # identify the step within each cycle and assign them a phase
  cycleorder <-
    cycleorder %>%
    filter((peak != 0) | (zerocross != 0)) %>%
    group_by(bodyparts, .add=TRUE) %>%
    mutate(cycle_step = case_when(zerocross == 1  ~  'up',
                                  peak == 1  ~  'max',
                                  zerocross == -1  ~  'down',
                                  peak == -1  ~  'min'),
           cycle_step = factor(cycle_step, levels = c('up', 'max', 'down', 'min')),
           tph = case_when(peak != 0  ~  tp,
                           zerocross != 0  ~  t0),
           ph1 = case_when(cycle_step == 'up'  ~  0,
                           cycle_step == 'max'  ~  0.25,
                           cycle_step == 'down'  ~  0.5,
                           cycle_step == 'min'  ~  0.75)) 
  
  if (!include.zeros) {
    cycleorder <-
      cycleorder %>%
      filter(cycle_step %in% c('min', 'max'))
  }
  
  # identify cycles. A cycle is defined as any time the phase steps backward. Each time it does
  # that, we increment the cycle number.
  cycleorder <-
    cycleorder %>%
    arrange(bodyparts, frame) %>%
    group_by(bodyparts) %>%
    mutate(cycle_num_step = as.integer(as.integer(cycle_step) < as.integer(lag(cycle_step))),
           cycle_num_step = replace_na(cycle_num_step, 0),
           cycle_num = cumsum(cycle_num_step),
           ph1 = ph1 + cycle_num) %>%
    select(bodyparts, frame, cycle_step, tph, cycle_num, ph1, peaksize, peakgap)

  # check and make sure we have at least two defined phases
  if (sum(!is.na(cycleorder$ph1)) < 2) {
    df <-
      df %>%
      ungroup() %>%
      mutate(freq = NA_real_,
             phase = NA_real_,
             cycle = NA_real_)

      trial1 <- df %>%
        distinct(trial) %>%
        slice_head(n = 1) %>%
        pull(trial)

      warning(sprintf('Cannot estimate frequency in trial %s', trial1))
  }
  else {
    if (all(is.na(filter(cycleorder, bodyparts == 'tailtip')$tph))) {
      tail = 'tailbase'
    } else {
      tail = 'tailtip'
    }
    
    # using a natural smoothing spline to estimate the phase at any time
    sf_tail <- with(filter(cycleorder, bodyparts == tail),
                    approxfun(tph, ph1, yleft = NA_real_, yright = NA_real_))

    # estimate phase and frequency (frequency is the first derivative of phase)
    df <-
      df %>%
      ungroup() %>%
      left_join(cycleorder, by = c("bodyparts", "frame")) %>%
      arrange(bodyparts, t) %>%
      group_by(bodyparts) %>%
      mutate(phase = sf_tail(t),
             freq = (lead(phase) - lag(phase)) / (lead(t) - lag(t)),
             cycle = floor(phase*2) / 2) %>%
      ungroup()
    
  }
  
  df
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

interpolate_amplitude_even <- function(df, s.even) {
  if (nrow(df) != length(s.even)) {
    df$s.even = NA
    df$amp.even = NA
  }
  else {
    # s.even[length(s.even)] <- min(s.even[length(s.even)], max(df$s))
    
    amp.even <-
      with(df, {
        good <- !is.na(s) & !is.na(amp)
        if ((sum(good) > 3) & good[1] & good[length(good)]) {
          spline(s[good], amp[good], xout = s.even)$y
        } else {
          NA
        }
    })
    
    df$amp.even = amp.even
    df$s.even = s.even
  }
  
  df
}

get_amplitude_envelope <- function(df) {
  s.even <- pracma::linspace(0, first(df$len), n = length(levels(df$bodyparts)))
  
  df %>%
    ungroup() %>%
    filter(!is.na(amp)) %>%
    arrange(trial, cycle, bodyparts) %>%
    group_by(trial, cycle) %>%
    group_modify(~ interpolate_amplitude_even(.x, s.even)) %>%
    mutate(swimvel = mean(swimvel, na.rm = TRUE)) %>%
    ungroup() %>%
    select(Treatment, trial, cycle, bodyparts, s.even, amp.even, swimvel)
}

get_all_kinematics <- function(df,
                               widthdata,
                               smooth.excursion = 0.2,
                               min.peak.gap = 0.01,
                               min.peak.size = 0.01) {
  df %>%
    distinct(trial) %>%
    slice_head(n = 1) %>%
    pull(trial) %>%
    print()
  
  df %>%
    get_arc_length() %>%
    interpolate_width(widthdata) %>%
    get_center_of_mass() %>%
    get_swim_vel_dir() %>%
    get_body_axis() %>%
    get_excursions() %>%
    get_curvature() %>%
    get_cycles(smooth.excursion = smooth.excursion,
               min.peak.gap = min.peak.gap,
               min.peak.size = min.peak.size) %>%
    get_amplitudes() %>%
    get_wavespeed() %>%
    get_wavelength()
}



