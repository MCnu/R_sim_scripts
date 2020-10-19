

kurzplot <- function (ID, SILENCE = F, linearized = F, 
                      msd_input_dir = paste(pull_dir,"data/ANALYTICAL_ENVIRONMENT/CALCULATED_MSD/",sep = ""),
                      merge_output_dir = paste(pull_dir,"data/ANALYTICAL_ENVIRONMENT/MSD_merged/",sep = "")) {
  library(tidyverse)
  library(stats)
  library(data.table)
  library(ggplot2)
  setwd(msd_input_dir)
  filenames <- list.files(pattern = ID)
  len <- length(filenames)
  DF <-
    fread(
      paste(
        msd_input_dir,
        filenames[1],
        sep = "/"
      )
    )
  n <- 2
  while (n < len + 1) {
    df <-
      fread(
        paste(
          msd_input_dir,
          filenames[n],
          sep = "/"
        )
      )
    #print(filenames[n])
    DF <-
      rbind(DF[, c("MSD", "taus", "experiment")], df[, c("MSD", "taus", "experiment")])
    n <- n + 1
  }
  mmModel <-
    nls(data = DF,
        MSD ~ L * taus ^ a,
        df,
        start = list(L = 0.01, a = 0.45))
  sumModel <- summary(mmModel)
  mod <- coef(mmModel)
  summary(mmModel, correlation = TRUE)
  G <- round(mod[1], 3)
  A <- round(mod[2], 3)
  limX <- max(DF$taus) / 3
  aves <- aggregate(DF[, 1], list(DF$taus), mean)
  colnames(aves) <- c("taus", "aveMSD")
  aves$predict <- mod[1] * aves$taus ^ mod[2]
  R <- (round(cor(aves$aveMSD, aves$predict), 6)) ^ 2
  SDs <- aggregate(DF[, 1], list(DF$taus), sd)
  colnames(SDs) <- c("taus", "SD")
  aves$SDadd <- aves$aveMSD + (SDs$SD)
  aves$SDsub <- aves$aveMSD - (SDs$SD)
  
  setwd(merge_output_dir)
  write.csv(
    DF,
    file = paste(ID, "MSDs merged.csv"),
    row.names = TRUE,
    quote = FALSE
  )
  merged <- fread(paste(ID, "MSDs merged.csv"))
  
  if (SILENCE == F) {
    print(paste("Gamma:", G, sep = " "))
    print(paste("Alpha:", A, sep = " "))
    print(paste("R squared:", R, sep = " "))
    
    cplot <- ggplot(aves, aes(x = log10(taus), y = log10(aveMSD))) +
      geom_smooth(
        data = merged,
        aes(x = log10(taus), y = log10(MSD)),
        col = "blue",
        fill = "blue",
        alpha = 0.7,
        size = 0.7
      ) +
      
      geom_point(data = DF,
                 aes(x = log10(taus), y = log10(MSD)),
                 alpha = 0.08) +
      geom_line(data = aves, aes(x = log10(taus), y = log10(predict)), col = "red") +
      theme(legend.position = "none") +
      theme(
        text = element_text(
          size = 18,
          colour = "black"
        ),
        axis.text = element_text(size = 16, color = "black")
      ) +
      coord_cartesian(ylim = c(-3,-0.5), xlim = c(-0.65, 1.35)) +
      labs(x = "Log10 Tau (s)",
           y = "Log10 MSD" ~ (µm ^ 2),
           subtitle = ID) +
      #annotate("text", x = 4, y= 0.08, label = paste("italic(R) ^ 2 ==", R), parse = TRUE, size = 10)+
      #annotate("text", x = 4, y = 0.09, label = paste("MSD = ", G , "x Tau^", A), size = 10)+
      theme()
    
    #png(paste(ID, ".png", sep=""), height = 1500, width = 1500, res = 300)
    print(cplot)
    
    
    #print(SDs)
    stud <<- SDs
    astu <<- aves
  }
  
  if(linearized == T){
    return(c(G,A,R))
  }
  
}