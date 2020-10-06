library(stats)
# library(boot)
library(data.table)
# library(DescTools)
library(tidyverse)
# library(mixtools)
library(viridis)

uap = "U3alonePlus"
gup = "G1atU3Plus"

comtheme <- theme_bw()+
            theme(text = element_text(color = "#545454", size = 18),
                  legend.position = "none")

#go to the analytical working directory
setwd("C:/USERS/MCS/Desktop/PULLS/R_sim_scripts/data/Analytical_environment")
list.files()




#form list of specific raw data as matricies, perform MSD analysis, plot,
#and merge MSDs for condition/sim
setwd("C:/USERS/MCS/Desktop/PULLS/R_sim_scripts/data/Analytical_environment/raw_data")
list.files()
setwd(list.files()[3])

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

simMSD(msd_list, "PY_URA3_20_70_FLE")
kurzplot("PY_URA3_20_70_FLE")

mmsdplot(c(uap, gup, "PY_URA3_00_00_FLE", "PY_URA3_20_70_FLE","PY_URA3_20_80_FLE", "PY_URA3_80_90_FLE","PY_URA3_90_70_FLE", "PY_URA3_90_75_FLE","PY_URA3_90_80_FLE","PY_SPB_100_100_FLE"))

#cut length of sim_list to remove buffer (300 secs or 1500 steps)

debuff_sim_list <- vector("list", length = length(sim_list))
debuff_end <- nrow(sim_list[[1]])
debuff_start <- debuff_end - 3000

for(i in 1:length(sim_list)){
  debuff_sim_list[[i]] <- sim_list[[i]][debuff_start:debuff_end,]
  
}

start_x <- c()
start_y <- c()
for(i in 1:length(debuff_sim_list)){
  start_x <- c(start_x,debuff_sim_list[[i]][1,2])
  start_y <- c(start_y,debuff_sim_list[[i]][1,1])
}

data.frame(x = start_x, y = start_y)

plot(x=start_x, y= start_y)


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
tot_d2o <- c()
d2_list <- vector("list", length(TBF))

for (i in 1:length(TBF)){
  float_mat <- TBF[[i]]
  d2ovec[1] <- sqrt(float_mat[1,1]^2 + float_mat[1,2]^2)
  for(j in 2:nrow(float_mat)){
    d2svec[j] <- sqrt((float_mat[1,1] - float_mat[j,1])^2 + (float_mat[1,2] - float_mat[j,2])^2)
    d2pvec[j] <- sqrt((float_mat[(j-1),1] - float_mat[j,1])^2 + (float_mat[(j-1),2] - float_mat[j,2])^2)
    d2ovec[j] <- sqrt(float_mat[j,1]^2 + float_mat[j,2]^2)
  }
  float_mat <-cbind(float_mat, matrix(d2svec, ncol = 1), matrix(d2pvec, ncol = 1), matrix(d2ovec, ncol =1))
  colnames(float_mat) <- c("x","y","is_bound","time", "D2P", "D2S", "D2O")
  row.names(float_mat) <- c(1:nrow(float_mat))
  d2_list[[i]] <- float_mat
  tot_d2o <- c(tot_d2o,d2ovec)
}

summary(tot_d2o)

plot(density(tot_d2o, from = 0, to = 1))

#Repositioning
#REPOSITIONING USES A NON_DEBUFF SIM_LIST passed through the d2_list

Repos_List <- d2_list
sim_lin()

tvv_frame <- rbind(tvv_frame, data.frame(velocity <- tvv, ID = rep("u9b8",length(tvv))))
mdt_frame <- rbind(mdt_frame, data.frame(mean_delta_theta = global_summary_sim[,"mean_delta_theta"], ID = rep("u9b8",nrow(global_summary_sim))))
arrival_frame <- rbind(arrival_frame, data.frame(arrival_time = global_summary_sim[,"time"], ID = rep("u9b8",nrow(global_summary_sim))))

ggplot(data = tvv_frame)+
  geom_violin(aes(x = tvv_frame[,2], y = tvv_frame[,1], fill = tvv_frame[,2]), size = 2, draw_quantiles = c(0.25, 0.5,0.75))+
  scale_fill_viridis(discrete = T)+
  comtheme+
  xlab("Simulation")+
  ylab("velocity (um/s)")

ks.test(tvv_frame[tvv_frame[,2]=="u2b8",1],tvv_frame[tvv_frame[,2]=="u9b8",1] )

ggplot(data = mdt_frame)+
  geom_violin(aes(x = mdt_frame[,2], y = mdt_frame[,1], fill = mdt_frame[,2]), size = 2, draw_quantiles = c(0.25, 0.5,0.75))+
  scale_fill_viridis(discrete = T)+
  comtheme+
  xlab("Simulation")+
  ylab("mean delta theta")


ks.test(mdt_frame[mdt_frame[,2]=="u2b8",1],mdt_frame[mdt_frame[,2]=="u9b8",1] )

ggplot(data = arrival_frame)+
  geom_violin(aes(x = arrival_frame[,2], y = arrival_frame[,1], fill = arrival_frame[,2]), size = 2, draw_quantiles = c(0.25, 0.5,0.75))+
  scale_fill_viridis(discrete = T)+
  comtheme+
  xlab("Simulation")+
  ylab("Time Before First Contact")

ks.test(arrival_frame[arrival_frame[,2]=="u2b8",1],arrival_frame[arrival_frame[,2]=="u9b8",1] )

#Periphery Occupancy

periph_occu_list <- vector("list", length(d2_list)) 
bound_zone_thickness <- (1-sqrt(2/3)) # this is outer third of nucleus
for(i in 1:length(d2_list)){
  float_po_logic <- d2_list[[i]][,"D2O"] > (1-bound_zone_thickness)
  float_po_matrix <- matrix(seq(0,(length(float_po_logic)-1)*.21,.21), ncol = 1)
  periph_occu_list[[i]] <- cbind(float_po_matrix, matrix(as.numeric(float_po_logic), ncol = 1))
}

preproc_NAs <- rep(NA,(length(periph_occu_list[[i]][,1])*length(periph_occu_list)))
po_plot_frame <- data.frame(Trajectory = preproc_NAs, time = preproc_NAs, is_periph = preproc_NAs)
for(i in 1:length(periph_occu_list)){
  po_plot_frame[((i-1) * nrow(periph_occu_list[[i]]) + 1):(i*nrow(periph_occu_list[[i]])),] <- data.frame(Trajectory = i, time = periph_occu_list[[i]][,1], is_periph = periph_occu_list[[i]][,2])
}

#po_plot <- ggplot(data = po_plot_frame)+
#  geom_tile(aes(x = time, y = Trajectory, color = is_periph), position = "dodge")+
#  scale_colour_viridis()

#print(po_plot)

size_of_bin = 50

preproc_NAs <- rep(NA,(length(po_plot_frame[,1])/size_of_bin))
bin_po_plot_frame <- data.frame(Trajectory = preproc_NAs, time = preproc_NAs, is_periph = preproc_NAs)
for(i in 1:length(bin_po_plot_frame[,1])){
  bin_po_plot_frame[i,1] <- mean(po_plot_frame[((i-1)*size_of_bin+1):(i*size_of_bin), 1])
  bin_po_plot_frame[i,2] <- po_plot_frame[((i-1)*size_of_bin+1), 2]
  bin_po_plot_frame[i,3] <- round(mean(po_plot_frame[((i-1)*size_of_bin+1):(i*size_of_bin), 3]))
}

#bpo_plot <- ggplot(data = bin_po_plot_frame)+
#  geom_tile(aes(x = time, y = Trajectory, color = is_periph), position = "dodge")+
#  scale_colour_viridis()

#print(bpo_plot)

bpo_boundfrac_traj <- c()
for(i in 1:length(unique(bin_po_plot_frame[,1]))){
  bpo_boundfrac_traj[i] <- mean(bin_po_plot_frame[(bin_po_plot_frame[,1] == unique(bin_po_plot_frame[,1])[i]),3])
}


#plot(density(bpo_boundfrac_traj))

#plot(y=bpo_boundfrac_traj,x=c(1:length(bpo_boundfrac_traj)))

#get fancy, reorder the bin_po_plot_frame by is_bound traj percentage

ordered_bpopf <- bin_po_plot_frame
ordered_bpopf <- mutate(ordered_bpopf, rank = NA)
for(i in 1:length(bpo_boundfrac_traj)){
  float_frame <- filter(bin_po_plot_frame, Trajectory == (order(bpo_boundfrac_traj)[i]))
  ordered_bpopf[((i-1)*nrow(float_frame) + 1):(i*nrow(float_frame)),1:3] <- float_frame
  ordered_bpopf[((i-1)*nrow(float_frame) + 1):(i*nrow(float_frame)),4] <- i
}

obpo_plot <- ggplot(data = ordered_bpopf)+
  geom_tile(aes(x = time, y = rank, color = is_periph, fill = is_periph), position = "dodge")+
  theme(legend.position = "none")+
  scale_colour_viridis()+
  scale_fill_viridis()+
  comtheme

print(obpo_plot)

bpo_boundfrac_time <- c()
for(i in 1:length(unique(bin_po_plot_frame[,2]))){
  bpo_boundfrac_time[i] <- mean(bin_po_plot_frame[(bin_po_plot_frame[,2] == unique(bin_po_plot_frame[,2])[i]),3])
}

plot(density(bpo_boundfrac_time))

plot(y=bpo_boundfrac_time,x=c(1:length(bpo_boundfrac_time)))

bpobf_time_frame <- data.frame(bound_frac = bpo_boundfrac_time, time = c(1:length(bpo_boundfrac_time)))

ggplot(data = bpobf_time_frame)+
  stat_smooth(aes(x= time, y = bound_frac*100))+
  geom_point(aes(x= time, y = bound_frac*100), alpha = 0.5)+
  coord_cartesian(ylim = c(0,100))+
  comtheme+
  ylab("% of Trajectories at the periphery")+
  xlab("time (s)")


#TODO:create comparative plot for u2b8 and u9b8 FLE sims

#u2b8_bpo_bf_traj <- bpo_boundfrac_traj
#u2b8_bpo_bf_time <- bpo_boundfrac_time

#u9b8_bpo_bf_traj <- bpo_boundfrac_traj
#u9b8_bpo_bf_time <- bpo_boundfrac_time

#u9b9_bpo_bf_traj <- bpo_boundfrac_traj
#u9b9_bpo_bf_time <- bpo_boundfrac_time
# 
# u0b0_bpo_bf_traj <- bpo_boundfrac_traj
# u0b0_bpo_bf_time <- bpo_boundfrac_time

comp_traj_frame <- data.frame(trajperc = c(u0b0_bpo_bf_traj,u2b8_bpo_bf_traj,u9b8_bpo_bf_traj,u9b9_bpo_bf_traj), 
                              ID = c(rep("u00/b00", length(u0b0_bpo_bf_traj)),
                                     rep("u20/b80", length(u2b8_bpo_bf_traj)), 
                                     rep("u90/b80", length(u9b8_bpo_bf_traj)),
                                     rep("u90/b90", length(u9b9_bpo_bf_traj))))


ctjplot <- ggplot(data=comp_traj_frame)+
  geom_violin(aes(x = ID, y = (trajperc*100), fill = ID), draw_quantiles = c(0.25, 0.5, 0.75), width = 5, size = 1, position = position_dodge(width = 20))+
  geom_point(aes(x=ID,y=(trajperc*100)), alpha = 0.2, position = position_jitter(height = 0, width = 0.1))+
  scale_fill_viridis(discrete = T)+
  scale_color_viridis(discrete = T)+
  coord_cartesian(ylim = c(0,100))+
  comtheme+
  ylab("% time at the periphery per Trajectory")+
  xlab("Unbound to Bound rate / Bound to Bound rate")

print(ctjplot)

comp_time_frame <- data.frame(timeperc = c(u0b0_bpo_bf_time,u2b8_bpo_bf_time,u9b8_bpo_bf_time,u9b9_bpo_bf_time), 
                              ID = c(rep("u00/b00", length(u0b0_bpo_bf_time)),
                                     rep("u20/b80", length(u2b8_bpo_bf_time)), 
                                     rep("u90/b80", length(u9b8_bpo_bf_time)),
                                     rep("u90/b90", length(u9b9_bpo_bf_time))))

ctmplot <- ggplot(data=comp_time_frame)+
  geom_violin(aes(x = ID, y = timeperc*100, fill = ID),size = 1, width = 3, draw_quantiles = c(0.25,0.5,0.75))+
  geom_point(aes(x=ID,y=(timeperc*100)), alpha = 0.05, position = position_jitter(height = 0, width = 0.1))+
  scale_fill_viridis(discrete = T)+
  coord_cartesian(ylim = c(0,100))+
  comtheme+
  ylab("% of Trajectories at periphery per second")

print(ctmplot)

#cluster trajectory comparison using d2_list from above

#TODO:turn d2_list into two cluster pair lists

d2_list




