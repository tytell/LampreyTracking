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