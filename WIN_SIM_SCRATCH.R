library(stats)
# library(boot)
library(data.table)
# library(DescTools)
library(tidyverse)
# library(mixtools)
library(viridis)

#go to the analytical working directory
setwd("C:/USERS/MCS/Desktop/Analytical_environment")
list.files()




#form list of specific raw data as matricies, perform MSD analysis, plot,
#and merge MSDs for condition/sim
setwd("C:/USERS/MCS/Desktop/Analytical_environment/raw_data")
list.files()
setwd(list.files()[6])

for(i in 1:length(list.files())){
  if(i == 1){
    sim_list <- vector(mode = "list", length = length(list.files()))
  }
  sim_list[[i]] <- as.matrix(fread(list.files()[i]))
}

for(i in 1:length(sim_list)){
  if(i == 1){
    msd_list <- vector(mode = "list", length = length(sim_list))
  }
  sl_len <- length(sim_list[[i]][,1])
  msd_list[[i]] <- sim_list[[i]][((sl_len - 199):sl_len),] 
}


tail(msd_list[[1]])

simMSD(msd_list, "PY_URA3_98_FLE")
kurzplot("PY_URA3_98_FLE")

mmsdplot(c(uap, gup, "URA3_29_FLE", "PY_URA3_28_FLE", "PY_URA3_98_FLE"))


 #Add time and "distance to" columns (start,previous,origin(aka radius))

TBF <- sim_list
subTBF_nrow <- length(TBF[[1]][,1])
time_seq <- matrix(seq(0, (subTBF_nrow-1)*0.21, 0.21), ncol = 1)
for(i in 1:length(TBF)){
  TBF[[i]] <- cbind(TBF[[i]],time_seq)
  colnames(TBF[[i]]) <- c("x","y","is_bound","time")
}


t_int <- TBF[[1]][,"time"][2]
d2svec <- c(0)
d2pvec <- c(0)
d2ovec <- c(0)
d2_list <- vector("list", length(TBF))

for(i in 1:length(TBF)){
  float_mat <- TBF[[i]]
  for(j in 2:nrow(float_mat)){
    d2svec[j] <- sqrt((float_mat[1,1] - float_mat[j,1])^2 + (float_mat[1,2] - float_mat[j,2])^2)
    d2pvec[j] <- sqrt((float_mat[(j-1),1] - float_mat[j,1])^2 + (float_mat[(j-1),2] - float_mat[j,2])^2)
    d2ovec[j] <- sqrt(float_mat[j,1]^2 + float_mat[j,2]^2)
  }
  float_mat <-cbind(float_mat, matrix(d2svec, ncol = 1), matrix(d2pvec, ncol = 1), matrix(d2ovec, ncol =1))
  colnames(float_mat) <- c("x","y","is_bound","time", "D2P", "D2S", "D2O")
  row.names(float_mat) <- c(1:nrow(float_mat))
  d2_list[[i]] <- float_mat
}

#Repositioning

Repos_List <- d2_list
sim_lin()


#Periphery Occupancy

periph_occu_list <- vector("list", length(d2_list)) 
bound_zone_thickness <- 0.15
for(i in 1:length(d2_list)){
  float_po_logic <- d2_list[[i]][,"D2O"] > (1-bound_zone_thickness)
  float_po_matrix <- matrix(seq(0,(length(float_po_logic)-1)*.21,.21), ncol = 1)
  periph_occu_list[[i]] <- cbind(float_po_matrix, matrix(as.numeric(float_po_logic), ncol = 1))
}

preproc_NAs <- rep(NA,(length(periph_occu_list[[i]][,1])*length(periph_occu_list)))
po_plot_frame <- data.frame(Trajectory = preproc_NAs, time = preproc_NAs, is_periph = preproc_NAs)
for(i in 1:length(periph_occu_list)){
  po_plot_frame[((i-1) * 4500 + 1):(i*4500),] <- data.frame(Trajectory = i, time = periph_occu_list[[i]][,1], is_periph = periph_occu_list[[i]][,2])
}

po_plot <- ggplot(data = po_plot_frame)+
  geom_tile(aes(x = time, y = Trajectory, color = is_periph), position = "dodge")+
  scale_colour_viridis()

print(po_plot)

preproc_NAs <- rep(NA,(length(po_plot_frame[,1])/5))
bin_po_plot_frame <- data.frame(Trajectory = preproc_NAs, time = preproc_NAs, is_periph = preproc_NAs)
for(i in 1:length(bin_po_plot_frame[,1])){
  bin_po_plot_frame[i,1] <- mean(po_plot_frame[((i-1)*5+1):(i*5), 1])
  bin_po_plot_frame[i,2] <- po_plot_frame[((i-1)*5+1), 2]
  bin_po_plot_frame[i,3] <- round(mean(po_plot_frame[((i-1)*5+1):(i*5), 3]))
}

bpo_plot <- ggplot(data = bin_po_plot_frame)+
  geom_tile(aes(x = time, y = Trajectory, color = is_periph), position = "dodge")+
  scale_colour_viridis()

print(bpo_plot)

po_boundfrac_traj <- c()
for(i in 1:length(unique(bin_po_plot_frame[,1]))){
  po_boundfrac_traj[i] <- mean(bin_po_plot_frame[(bin_po_plot_frame[,1] == unique(bin_po_plot_frame[,1])[i]),3])
}

po_boundfrac_time <- c()
for(i in 1:length(unique(bin_po_plot_frame[,2]))){
  po_boundfrac_time[i] <- mean(bin_po_plot_frame[(bin_po_plot_frame[,2] == unique(bin_po_plot_frame[,2])[i]),3])
}
po_boundfrac_time

plot(density(po_boundfrac_traj))
plot(density(po_boundfrac_time))

