

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

sim_clust_magcor <- function(primary_xy, secondary_xy) {
  plen <- length(primary_xy)
  slen <- length(secondary_xy)
  if (plen != slen) stop("lists must have equal length")

  sub_plen <- nrow(primary_xy[[1]])
  sub_slen <- nrow(secondary_xy[[1]])
  if (sub_plen != sub_slen) stop("unequal d2p sub matrix")

  cor_output_mat <- matrix(data = NA, nrow = plen, ncol = 2)

  for (i in 1:plen) {
    # initiate d2p matrix
    d2p_matrix <- matrix(data = NA, nrow = sub_plen - 1, ncol = 2)

    for (j in 2:sub_plen) {
      pri_x <- (primary_xy[[i]][j, 1] - primary_xy[[i]][(j - 1), 1])
      sec_x <- (secondary_xy[[i]][j, 1] - secondary_xy[[i]][(j - 1), 1])
      pri_y <- (primary_xy[[i]][j, 2] - primary_xy[[i]][(j - 1), 2])
      sec_y <- (secondary_xy[[i]][j, 2] - secondary_xy[[i]][(j - 1), 2])
      p_leg <- sqrt(pri_x^2 + pri_y^2)
      s_leg <- sqrt(sec_x^2 + sec_y^2)
      d2p_matrix[j - 1, ] <- c(p_leg, s_leg)
    }
    float_cor <- cor(d2p_matrix[, 1], d2p_matrix[, 2])
    cor_output_mat[i, c(1, 2)] <- c(i, float_cor)
  }

  cor_output_df <- as.data.frame(cor_output_mat)
  colnames(cor_output_df) <- c("Pair_ID", "mag_PCC")

  return(cor_output_df)
}

sim_clust_theta <- function(primary_xy, secondary_xy) {
  plen <- length(primary_xy)
  slen <- length(secondary_xy)
  if (plen != slen) stop("lists must have equal length")

  sub_plen <- nrow(primary_xy[[1]])
  sub_slen <- nrow(secondary_xy[[1]])

  if (sub_plen != sub_slen) stop("unequal d2p sub matrix")

  theta_output_list <- lapply(1:plen, matrix, data = NA, nrow = sub_plen, ncol = 1)

  for (i in 1:plen) {
    for (j in 2:sub_plen) {
      pri_x <- (primary_xy[[i]][j, 1] - primary_xy[[i]][(j - 1), 1])
      sec_x <- (secondary_xy[[i]][j, 1] - secondary_xy[[i]][(j - 1), 1])
      pri_y <- (primary_xy[[i]][j, 2] - primary_xy[[i]][(j - 1), 2])
      sec_y <- (secondary_xy[[i]][j, 2] - secondary_xy[[i]][(j - 1), 2])
      p_leg <- sqrt(pri_x^2 + pri_y^2)
      s_leg <- sqrt(sec_x^2 + sec_y^2)
      i_leg <- sqrt((pri_x - sec_x)^2 + (pri_y - sec_y)^2)

      if (p_leg < 0.001 | s_leg < 0.001) {
        next
      }

      float_acos_arg <- (i_leg^2 - p_leg^2 - s_leg^2) / (-1 * 2 * p_leg * s_leg)

      float_ang <- acos(float_acos_arg)

      if (is.nan(float_ang)) {
        print(c(float_acos_arg, p_leg, s_leg, i_leg))
        next
      } else {
        theta_output_list[[i]][(j - 1), 1] <- float_ang
      }
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
