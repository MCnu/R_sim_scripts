comtheme <- theme(text = element_text(family = "Arial", color = "#545454", size = 21),
                  axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  legend.position = "none")


mmsdplot <- function(IDarray, spbmsd = T, title="test", plottitle= substitute(paste(italic("URA3"))), merge_input_dir = "C:/USERS/MCS/DESKTOP/ANALYTICAL_ENVIRONMENT/MSD_merged", plot_output_dir, save_plot = F){
  
  
  setwd(merge_input_dir)
  msdcomplete <- data.frame(MSD = numeric(), taus = numeric(), experiment = character(), ID = character())
  for(i in 1:length(IDarray)){
    fname <- list.files(pattern = IDarray[i])
    print(fname)
    dffloat <- fread(list.files(pattern = IDarray[i]))
    #print(dffloat[,c("MSD","taus","experiment")])
    msdcomplete <- rbind(msdcomplete, cbind(dffloat[,c("MSD","taus","experiment")], data.frame(ID = as.character(IDarray[i]))))
    
  }
  
  print(msdcomplete)
  
  #saves plot axis and internal annotation label names
  xaxlabel = expression("Tau (us)")
  yaxlabel = expression("Mean squared displacement" ~ (µm^2))
  
    spb_float_x <- c(1:200)*0.21
    spb_float_y <- (spb_float_x * 0.0009) + 0.0015
    
    mplot <- ggplot()+
      comtheme+
      #Set up  log scale
      scale_x_log10(limits = c(0.2, 2.1))+
      scale_y_log10(limits = c(0.000999,0.1))+
      annotation_logticks(color = "black")+
      geom_smooth(data = msdcomplete, aes(x=taus, y=MSD, color = ID, fill = ID))+
      labs(x=xaxlabel, y = yaxlabel, title = plottitle)+
      #scale_fill_manual(values = c("#66A53D", "#52307c", "#BC6C45"))+
      #scale_color_manual(values = c("#66A53D", "#52307c", "#BC6C45"))+
      theme(legend.position = "right", plot.title = element_blank())+
      coord_cartesian()
    
    if(length(IDarray) == 2){
      mplot<- mplot+
        scale_fill_manual(values = c("#52307c", "#66A53D"))+
        scale_color_manual(values = c("#52307c", "#66A53D"))
    }else{
      mplot <- mplot + 
        scale_color_viridis(end = 0.8, discrete = T)+
        scale_fill_viridis(end = 0.8, discrete = T)
    }
    if(spbmsd == T){
      mplot <- mplot+
        geom_line(data=data.frame(xs = spb_float_x, ys = spb_float_y), aes(x=xs,y=ys), size = 2, alpha = 0.5, linetype = "dashed")
    }
    
    
  
 
  if(save_plot == T){
    #Save and Print the plot
    ggsave(plot_output_dir, plot = mplot, width = 6, height = 4.8)
  }
  
  print(mplot)
  print("Plot's done!")
  
  global_mmsd <<- msdcomplete
  
}
