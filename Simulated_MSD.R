

simMSD <- function(input_list, fileID, taus = 100, SILENCE = T,
                   output_dir = paste(dataset_dir, "trajecotry_r_data/analytical_environment/calculated_MSD/", sep = "")) {
  for (i in 1:length(input_list)) {
    DF <- as.data.frame(input_list[[i]])
    
    flen <- nrow(DF)
    time_interval <- DF$time[3] - DF$time[2]
  
    # REMAKE TO PROPERLY ESTABLISH DIMENSIONS
    MSDmat <- matrix(nrow = taus, ncol = 2)
    MSDmat[, 2] <- c(1:taus)
    for (j in 1:taus) {
      floatMSD <- matrix(nrow = (flen - j), ncol = 2)
      floatMSD[, 2] <- c(1:(flen - j))
      for (k in 1:(flen - j)) {
        MeSqDi <- ((DF[, 1][k] - DF[, 1][(k + j)])^2 + (DF[, 2][k] - DF[, 2][(k + j)])^2)
        floatMSD[k, 1] <- MeSqDi
      }
      MSDmat[j, 1] <- mean(floatMSD[, 1])
    }
    MSD <- data.frame(MSD = as.numeric(MSDmat[, 1]), taus = as.numeric(MSDmat[, 2] * time_interval), experiment = paste("SIM", i, sep = "_"))
    setwd(output_dir)
    write.csv(MSD, paste("MSD ", fileID, "_", i, ".csv", sep = ""), row.names = F)
    # print(MSD)
  }
  if (SILENCE == F) {
    print(paste("Saved ", length(input_list), " files for simulated paths.", sep = ""))
  }
}
