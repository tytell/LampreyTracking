require(tidyverse)

get_arc_length <- function(df) {
  df %>%
    group_by(frame) %>%
    mutate(dx = lead(xmm) - xmm,
           dy = lead(ymm) - ymm,
           s = cumsum(sqrt(dx^2 + dy^2)),
           s = replace_na(lag(s), 0)) %>%
    ungroup()
}

get_length <- function(df) {
  df %>%
    group_by(frame) %>%
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
  
  len <- get_length(df)
  
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
    group_by(frame) %>%
    summarize(m = pracma::trapz(s, pi*height*width/4),
              comx = pracma::trapz(s, xmm * pi*height*width/4) / m,
              comy = pracma::trapz(s, ymm * pi*height*width/4) / m)
  
  df %>%
    left_join(com, by = "frame")
}

smooth_com_spline <- function(x,y, spar) {
  ys <- numeric(length(y))
  good <- !is.na(y)
  
  sp <- smooth.spline(x[good], y[good], spar = spar)
  
  ys[!good] <- NA
  ys[good] = predict(sp)$y
  
  ys
}

get_swim_vel_dir <- function(df, s = 0.8, fps = 60) {
  swim <-
    df %>%
    group_by(frame) %>%
    summarize(comx = first(comx),
              comy = first(comy),
              t = first(t)) %>%
    mutate(swimvelx = (lead(comx) - lag(comx)) * fps / 2,
           swimvely = (lead(comy) - lag(comy)) * fps / 2,
           swimvel = sqrt(swimvelx^2 + swimvely^2),
           swimvelxs = smooth_com_spline(t, swimvelx, s),
           swimvelys = smooth_com_spline(t, swimvely, s),
           swimvels = sqrt(swimvelxs^2 + swimvelys^2),
           swimdirx = swimvelxs / swimvels,
           swimdiry = swimvelys / swimvels) %>%
    select(-comx, -comy, -t)
  
  df %>%
    left_join(swim, by = c("frame"))
}

get_excursions <- function(df) {
  df %>%
    mutate(exc = -(xmm - comx) * swimdiry + 
             (ymm - comy) * swimdirx)
}

get_cycles <- function(df) {
  df <-
    df %>%
    group_by(bodyparts) %>%
    mutate(excn = lead(exc),
           zerocross = case_when(sign(excn) > sign(exc)   ~   1,
                                 sign(excn) == sign(exc)  ~   0,
                                 sign(excn) < sign(exc)   ~   -1,
                                 TRUE   ~   0))
  
  zerocross <-
    df %>%
    group_by(bodyparts) %>%
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

  t0body <-
    df %>%
    ungroup() %>%
    filter(zerocross != 0 &
             bodyparts != "tailtip") %>%
    group_by(bodyparts) %>%
    mutate()
  df %>%
    ungroup() %>%
    mutate(freq = approx(t0tail$t0, t0tail$freqbody, t)$y,
           phase = approx(t0tail$t0, t0tail$phtail, t)$y,
           cycle = floor(phase*2) / 2)
}

get_amplitudes <- function(df) {
  amps <-
    df %>%
    group_by(bodyparts, cycle) %>%
    summarize(ampframe = frame[which.max(abs(exc))],
              amp = max(abs(exc))) %>%
    mutate(amp = (lag(amp) + 2*amp + lead(amp)) / 4) %>%
    ungroup()
  
  left_join(df, amps, by = c("bodyparts", "cycle", "frame" = "ampframe"))
}

get_next_level <- function(fct) {
  lev <- levels(fct)
  n = as.numeric(fct) + 1
  
  good = (n >= 1) & (n <= length(lev))
  n[!good] = 1
  
  nextfct <- factor(lev[n], levels = lev)
  nextfct[!good] = NA
  nextfct
}

get_next_t0 <- function(a,b, sgn, nextdf) {
  nextt0 <- 
    nextdf %>%
    filter((t0 >= a) & (t0 < b) & (zerocross == sgn)) %>%
    transmute(t0nextseg = t0,
           snextseg = s)
  
  if (nrow(nextt0) == 0) {
    nextt0 = data.frame(t0nextseg = NA, snextseg = NA)
  }
  nextt0[1,]
}

get_wavespeed <- function(df) {
  bodyparts <- levels(df$bodyparts)
  bodyparts <- factor(bodyparts, levels = bodyparts)

  t0seg = list()
  for (i in seq(length.out = length(bodyparts)-1)) {
    seg1 <- bodyparts[i]
    
    seg2 = get_next_level(seg1)
    
    t0seg1 <- data1 %>%
      ungroup() %>%
      filter(bodyparts == seg1 & !is.na(t0)) %>%
      arrange(t0) %>%
      mutate(t0next = lead(t0)) %>%
      select(bodyparts, t0, t0next, zerocross)
    
    t0seg2 = data1 %>%
      ungroup() %>%
      filter(bodyparts == seg2 & !is.na(t0)) %>%
      select(bodyparts, s, t0, zerocross)
    
    t0seg1$t0nextseg <- NA
    t0seg1$snextseg <- NA
    for (r in seq(length.out = nrow(t0seg1))) {
        nextt0 <- 
          t0seg2 %>%
          filter((t0 >= t0seg1$t0[r]) & (t0 < t0seg1$t0next[r]) & (zerocross == t0seg1$zerocross[r])) %>%
          transmute(t0nextseg = t0,
                    snextseg = s) %>%
          filter(!is.na(t0nextseg))
      
      t0seg1$t0nextseg[r] <- pull(nextt0, t0nextseg) %>% first()
      t0seg1$snextseg[r] <- pull(nextt0, snextseg) %>% first()
    }
    
    t0seg1$nextseg = seg2
    t0seg[[i]] = t0seg1
  }
  
  t0seg <- bind_rows(t0seg)
  
  df %>%
    left_join(t0seg, by = c("bodyparts", "t0", "zerocross"))
}
