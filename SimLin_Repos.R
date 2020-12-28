# SIMULATED REPOSITIONING ANALYSIS

sim_lin <-
  function(delta_start = "follow",
           time_min = 25,
           time_max = 450,
           peri_rad = 0.85,
           velo = T,
           JDD = F,
           Direct = T,
           vacf = T,
           ctrl_summ = F,
           summarize = T,
           contact = "first") {
    # create LUT
    LUT_ttp_vec <- c()
    for (i in 1:length(Repos_List)) {
      float_mat <- Repos_List[[i]]

      # (nrow(float_mat[float_mat[,"D2O"] > peri_rad,]))
      float_ttp <- float_mat[float_mat[, "D2O"] > peri_rad, ]
      if (!is.null(nrow(float_ttp))) {
        if (nrow(float_ttp) >= 1) {
          LUT_ttp_vec[i] <- float_ttp[1, "time"]
        } else {
          LUT_ttp_vec[i] <- NA
        }
      } else {
        LUT_ttp_vec[i] <- NA
      }
    }

    LUT_ttpv_logic <-
      (LUT_ttp_vec > time_min & LUT_ttp_vec < time_max & !is.na(LUT_ttp_vec))
    print("Fraction of used trajectories:")
    print(mean(LUT_ttpv_logic))
    
    JDD_DF <-
      data.frame(
        cID = character(),
        cont_num = numeric(),
        position = character(),
        displacement = numeric()
      )
    DIRECT_DF <- data.frame()
    pridist_vec <- c()

    if(velo==T) tot_velo_vec <- c()
    if(summarize==T) summ_frame <- data.frame()
    
    for (i in 1:length(LUT_ttpv_logic)) {
      if (!LUT_ttpv_logic[i]) {
        next
      }
      float_mat <- Repos_List[[i]]
      float_mat <- float_mat[float_mat[, "time"] <= LUT_ttp_vec[i], ]
      step_count <- nrow(float_mat)
      # print(precon_float_direct_matrix)
      float_direct_frame <- data.frame()
      if (delta_start == "follow") {
        follow_float_frame <- data.frame()
        endp <- float_mat[step_count, ]
        for (s in 2:step_count) {
          prip <- float_mat[(s - 1), ]
          secp <- float_mat[s, ]
          pridist <-
            sqrt(((prip[1]) - (endp[1]))^2 + ((prip[2]) - (endp[2]))^2)
          secdist <-
            sqrt(((secp[1]) - (endp[1]))^2 + ((secp[2]) - (endp[2]))^2)
          interdist <-
            sqrt(((prip[1]) - (secp[1]))^2 + ((prip[2]) - (secp[2]))^2)
          delta_dist_to_end <- unlist(pridist - secdist)
          if (s == 2) {
            pridist_vec <- c(pridist_vec, as.numeric(pridist))
          }
          # positive means it used to be further away
          if (pridist != secdist) {
            # law of cosines, baby
            deviation_angle <-
              unlist(acos((pridist^2 + interdist^2 - secdist^2) / (2 * pridist *
                interdist)))
            # print(deviation_angle)
            if (is.nan(deviation_angle)) {
              deviation_angle <- 0
            }
          }
          if (pridist == secdist) {
            deviation_angle <- 0
          }


          follow_float_frame <-
            rbind(
              follow_float_frame,
              data.frame(
                cID = i,
                cont_num = 1,
                step = (s - 1),
                delta_dist = delta_dist_to_end,
                delta_theta = deviation_angle
              )
            )
        }

        float_direct_frame <-
          rbind(float_direct_frame, follow_float_frame)

        if (velo == T) {
          # vector of velocities derived from the D2P and time per frame
          precontact_steps <- as.vector(float_mat[2:nrow(float_mat), "D2P"])
          time_per_frame <- float_mat[2,"time"] - float_mat[1,"time"]
          
          tot_velo_vec <-
            c(tot_velo_vec, c(precontact_steps / time_per_frame))

        }

        if (summarize == T) {
          summ_frame <-
            rbind(
              summ_frame,
              data.frame(
                cID = i,
                cont_num = "First",
                mean_delta_dist = mean(float_direct_frame$delta_dist),
                mean_delta_theta = mean(float_direct_frame$delta_theta),
                time = step_count * (float_mat[2, "time"] - float_mat[1, "time"])
              )
            )
        }
      }

      DIRECT_DF <- rbind(DIRECT_DF, float_direct_frame)
    }
    if (JDD == T) {
      sum_filt <-
        data.frame(
          cID = c(),
          cont_num = c(),
          position = c(),
          med_dis = c()
        )
      for (i in 1:length(unique(JDD_DF$cID))) {
        filt_jdf <- filter(JDD_DF, cID == unique(JDD_DF$cID)[i])
        filt_pre <- filter(filt_jdf, position == "PRE")
        filt_pos <- filter(filt_jdf, position == "POS")
        prebind <-
          data.frame(
            cID = filt_pre[1, 1],
            cont_num = filt_pre[1, 2],
            position = filt_pre[1, 3],
            med_dis = median(filt_pre[, 4])
          )
        posbind <-
          data.frame(
            cID = filt_pos[1, 1],
            cont_num = filt_pos[1, 2],
            position = filt_pos[1, 3],
            med_dis = median(filt_pos[, 4])
          )
        sum_filt <- rbind(sum_filt, prebind, posbind)
      }

      jplot <- ggplot() +
        geom_violin(data = JDD_DF, aes(x = cID, y = displacement, color = position)) +
        geom_hline(yintercept = sum_u3_d2p[3]) +
        geom_point(data = sum_filt, aes(x = cID, y = med_dis, color = position))

      # print(jplot)
    }
    if (delta_start == "follow") {
      global_direct_sim <<- DIRECT_DF
      pdvec <<- pridist_vec
    }
    if (velo == T) {
      tvv <<- tot_velo_vec
      # print(plot(density(tot_velo_vec)))
    }

    if (summarize == T) {
      global_summary_sim <<- summ_frame
      sumplot <- ggplot() +
        coord_cartesian(xlim = c(-0.01, 0.01), ylim = c(0, 3.14)) +
        geom_text(
          data = summ_frame,
          aes(x = mean_delta_dist, y = mean_delta_theta, label = time)
        )

      # print(sumplot)
    }
  }
