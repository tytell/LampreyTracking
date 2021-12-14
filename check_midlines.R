require(tidyverse)

get_point_order <- function(df) {
  df %>%
    group_by(trial, frame) %>%
    group_modify(~ get_point_order_in_frame(.x))
}

get_point_order_in_frame <- function(df) {
  # get the distances between all of the points in this frame
  xy <- select(df, xmm, ymm)
  D <- dist(xy)
  D <- as.matrix(D)
  
  # step through each point and get the next point that's closest
  n <- dim(D)[1]
  ind <- rep(NA_integer_, n)
  for (i in seq_len(n-1)) {
    d1 <- D[(i+1):n,i]
    if (any(!is.na(d1))) {
      ind[i] = which.min(d1) + i
    }
  }
  
  df$nextpoint <- ind
  
  df
}

spline_outlier_replace <- function(df) {
  if (sum(!df$isoutlier) < 3) {
    df$xmm <- NA_real_
    df$ymm <- NA_real_
  } else {
    t1 <- with(df, t[!isoutlier])
  
    span <- with(df, (t >= t1[1]) & (t <= t1[length(t1)]))
    
    df <-
      df %>%
      mutate(xmm = spline(t1, xmm[!isoutlier], xout = t)$y,
             ymm = spline(t1, ymm[!isoutlier], xout = t)$y)
    
    df$xmm[!span] = NA_real_
    df$ymm[!span] = NA_real_
  }
  
  df
}

runmed_check <- function(x, k, ...) {
  if (sum(!is.na(x)) < k) {
    xmed <- rep(NA_real_, length(x))
  } else {
    xmed <- runmed(x, k, ...)
  }
}

detect_and_replace_outliers <- function(df, k = 31, outlierthreshold = 30, 
                                   likelihoodthreshold = 0.4) {
  #' Removes outliers based on a median filter.
  #' 
  #' After https://ocefpaf.github.io/python4oceanographers/blog/2015/03/16/outlier_detection/
  #' Takes a running median of the data, then looks at the difference
  #' between each value and the running median at that point. Gets a median value of
  #' the difference across the entire data set, then looks for points where
  #' the cartesian (x,y) distance between a point and its running median value
  #' is greater than some multiple of the overall median.
  #' 
  #' @param df Data frame containing xmm and ymm coordinates.
  #' @param k Length of the running median filter (default 31 points).
  #' @param outlierthreshold Threshold to identify an outlier (default 30 times the baseline difference)
  #' @param likelihoodthreshold Threshold to remove a point based on DeepLabCut likelihood (default 0.4)
  
  df <-
    df %>%
    arrange(bodyparts, t) %>%
    group_by(bodyparts) %>%
    mutate(xmm0 = xmm,
           ymm0 = ymm,
           xmm1 = if_else(likelihood < likelihoodthreshold, NA_real_, xmm),
           ymm1 = if_else(likelihood < likelihoodthreshold, NA_real_, ymm),
           xmed = runmed_check(xmm1, k, endrule = 'median', na.action = "na.omit"),
           ymed = runmed_check(ymm1, k, endrule = 'median', na.action = "na.omit"),
           dxmed = xmm1 - xmed,
           dymed = ymm1 - ymed,
           dist_to_med = sqrt(dxmed^2 + dymed^2))
  
  df <-
    df %>%
    ungroup() %>%
    mutate(#med_dist_to_med = median(dist_to_med, na.rm = TRUE),
      outlierdist = dist_to_med, # / med_dist_to_med,
      isoutlier = outlierdist > outlierthreshold | is.na(outlierdist) | is.na(xmm1)) %>%
    select(-c(xmm1, ymm1, dxmed, dymed, xmed, ymed, dist_to_med))
  
  df %>%
    group_by(bodyparts) %>%
    group_modify(~ spline_outlier_replace(.x))
}

get_longest_block <- function(df, n.elim = 3, jump.dist = 10,
                              checkbodyparts = c('snout', 'tailtip')) {
  # look for jumps of greater than jump.dist and remove
  # those frames. Repeat n.elim times
  df2 <- df %>%
    filter(bodyparts %in% checkbodyparts)
  
  for (i in seq_len(n.elim)) {
    df2 <-
      df2 %>%
      group_by(bodyparts) %>%
      mutate(lagx = lag(x),
             lagy = lag(y),
             lagx = if_else(is.na(lagx), x, lagx),
             lagy = if_else(is.na(lagy), y, lagy)) %>%
      mutate(d = sqrt((x - lagx)^2 + (y - lagy)^2),
             isjump = d > jump.dist) %>%
      ungroup() %>%
      group_by(frame) %>%
      filter(all(!isjump)) %>%
      ungroup()
  }
  
  blocks <-
    df2 %>%
    group_by(bodyparts) %>%
    mutate(d = sqrt((x - lag(x))^2 + (y - lag(y))^2),
           isjump = d > jump.dist) %>%
    ungroup() %>%
    select(frame, isjump) %>%
    mutate(dblock = if_else(isjump & !lag(isjump, default = TRUE), 1, 0),
           dblock = if_else(is.na(dblock), 1, dblock),
           blocknum = cumsum(dblock)) %>%
    group_by(blocknum) %>%
    summarize(blocklen = n(),
              startframe = first(frame),
              endframe = last(frame)) %>%
    arrange(desc(blocklen))
  
  f1 <- slice_head(blocks, n = 1) %>%
    pull(startframe)
  f2 <- slice_head(blocks, n = 1) %>%
    pull(endframe)
  
  goodframes <-
    df2 %>%
    ungroup() %>%
    filter((frame >= f1) & (frame <= f2)) %>%
    distinct(frame) %>%
    pull(frame)
 
  df %>%
    filter(frame %in% goodframes)
 }