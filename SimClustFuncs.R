

sim_clust_interdist <- function(primary_xy, secondary_xy) {
  plen <- length(primary_xy)
  slen <- length(secondary_xy)
  if (plen != slen) stop("lists must have equal length")

  interdist_output_matlist <- lapply(1:plen, matrix, data = NA, nrow = nrow(primary_xy[[1]]), ncol = 1)

  for (i in 1:plen) {
    for (j in 1:nrow(primary_xy[[i]])) {
      float_dist <- sqrt((primary_xy[[i]][j, 1] - secondary_xy[[i]][j, 1])^2 + (primary_xy[[i]][j, 2] - secondary_xy[[i]][j, 2])^2)
      # This way is easier to read and can be directly plugged into the theta
      interdist_output_matlist[[i]][j, 1] <- float_dist
    }
  }
  return(interdist_output_matlist)
}

sim_clust_magcor <- function(primary_d2p, secondary_d2p) {
  plen <- length(primary_d2p)
  slen <- length(secondary_d2p)
  if (plen != slen) stop("lists must have equal length")

  cor_output_mat <- matrix(data = NA, nrow = plen, ncol = 2)

  for (i in 1:plen) {
    float_cor <- cor(as.vector(primary_d2p[[i]]), as.vector(secondary_d2p[[i]]))
    cor_output_mat[i, c(1, 2)] <- c(i, float_cor)
  }

  cor_output_df <- as.data.frame(cor_output_mat)
  colnames(cor_output_df) <- c("Pair_ID", "mag_PCC")

  return(cor_output_df)
}

sim_clust_theta <- function(primary_d2p, secondary_d2p, stepwise_interdist) {
  plen <- length(primary_d2p)
  slen <- length(secondary_d2p)
  ilen <- length(stepwise_interdist)
  if (plen != slen) stop("lists must have equal length")
  if (plen != ilen) stop("lists must have equal length")

  sub_plen <- length(primary_d2p[[1]])
  sub_slen <- length(secondary_d2p[[1]])
  sub_ilen <- nrow(stepwise_interdist[[1]])
  if (sub_plen != sub_slen) stop("unequal d2p sub matrix")
  if ((sub_plen + 1) != sub_ilen) stop("interdist should be one row greater than d2p")

  theta_output_list <- lapply(1:plen, matrix, data = NA, nrow = sub_plen, ncol = 1)

  for (i in 1:plen) {
    for (j in 1:sub_plen) {
      p_leg <- as.numeric(primary_d2p[[i]][j])
      s_leg <- as.numeric(secondary_d2p[[i]][j])
      i_leg <- as.numeric(stepwise_interdist[[i]][j + 1])

      # law of cosines
      float_acos_arg <- (-(i_leg^2 - p_leg^2 - s_leg^2)) / (2 * p_leg * s_leg)

      if (float_acos_arg > 1) float_acos_arg <- 1

      if (float_acos_arg < -1) float_acos_arg <- -1

      float_ang <- acos(float_acos_arg)

      if (is.nan(float_ang)) print(float_acos_arg)

      theta_output_list[[i]][j, 1] <- float_ang
    }
  }

  return(theta_output_list)
}

sim_clust_lifefrac <- function(interdist_matlist, time_btwn_steps) {
  ilen <- length(interdist_matlist)

  life_output_mat <- matrix(data = NA, nrow = ilen, ncol = 2)
  frac_output_mat <- matrix(data = NA, nrow = ilen, ncol = 2)

  for (i in 1:ilen) {
    clust_points_mat <- interdist_matlist[[i]][(interdist_matlist[[i]][, 1] <= 0.55), ]

    if (length(clust_points_mat) == 0) {
      life_output_mat[i, ] <- c(0, i)
      frac_output_mat[i, ] <- c(0, i)
      next
    }

    float_life <- length(clust_points_mat) * time_interval / 60
    float_frac <- length(clust_points_mat) / nrow(interdist_matlist[[i]])


    life_output_mat[i, ] <- c(float_life, i)
    frac_output_mat[i, ] <- c(float_frac, i)
  }

  output_list <- list(life_output_mat, frac_output_mat)

  return(output_list)
}
