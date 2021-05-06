library(tidyverse)
library(data.table)
library(viridis)
library(boot)

steFunc <- function(x, i) {
  sd(x[i]) / sqrt(length(x))
}

mmsdplot <-
  function(IDarray,
           spbmsd = T,
           title = "test",
           plot_header = NA,
           merge_input_dir = "./",
           plot_output_file,
           plot_legend = F,
           save_plot = F) {
    initial_dir <- getwd()
    setwd(merge_input_dir)
    msdcomplete <-
      data.frame(
        MSD = numeric(),
        taus = numeric(),
        experiment = character(),
        ID = character()
      )
    msdaves <- data.frame(aveMSD = numeric(), 
                          boot_ste = numeric(), 
                          taus = numeric(), 
                          ID = character()
                          )
    for (i in 1:length(IDarray)) {
      fname <- list.files(pattern = IDarray[i])
      print(fname)
      dffloat <- fread(list.files(pattern = IDarray[i]))
      # print(dffloat[,c("MSD","taus","experiment")])
      complete_float <- cbind(
        dffloat[, c("MSD", "taus", "experiment")],
        data.frame(ID = as.character(IDarray[i]))
      )
      exp_length <- length(unique(complete_float$experiment))
      print(exp_length)
      aves_float <- cbind(
        data.frame(aveMSD = rep.int(0, nrow(complete_float)), boot_ste = rep.int(0, nrow(complete_float))),
        dffloat[c(1:(nrow(complete_float) / exp_length)), "taus"],
        data.frame(ID = IDarray[i])
      )

      for (j in 1:nrow(aves_float)) {
        MSD_float <- unlist(complete_float[complete_float$taus == aves_float[j, "taus"], "MSD"])
        aves_float[j, "aveMSD"] <- mean(MSD_float)
        boot_obj_float <- boot(MSD_float, steFunc, 100)
        aves_float[j, "boot_ste"] <- boot_obj_float$t0
      }


      msdcomplete <-
        rbind(msdcomplete, cbind(
          dffloat[, c("MSD", "taus", "experiment")],
          data.frame(ID = as.character(IDarray[i]))
        ))
      msdaves <- rbind(msdaves, aves_float)
      
    }
    
    # saves plot axis and internal annotation label names
    xaxlabel <- expression("Tau (s)")
    yaxlabel <- expression("Mean squared displacement" ~ (μm^2))

    spb_float_x <- c(1:200) * 0.21
    spb_float_y <- (spb_float_x * 0.0009) + 0.0015

    breaks <- 10^(-10:10)
    minor_breaks <- rep(1:9, 21) * (10^rep(-10:10, each = 9))


    mplot <- ggplot() +
      comtheme +
      # Set up  log scale
      scale_x_log10(breaks = breaks, minor_breaks = minor_breaks, limits = c(0.2, 4.2)) +
      scale_y_log10(breaks = breaks, minor_breaks = minor_breaks, limits = c(0.000999, 0.1)) +
      annotation_logticks(color = "black") +
      geom_path(data = msdaves, aes(
        x = taus,
        y = aveMSD,
        color = ID,
      )) +
      geom_ribbon(data = msdaves, aes(
        x = taus,
        ymin = aveMSD - boot_ste, 
        ymax = aveMSD + boot_ste, 
        color = ID,
        fill = ID
      ),alpha = 0.5, color = NA) +
      labs(x = xaxlabel, y = yaxlabel) +
      coord_cartesian()

    if (class(plot_header) == "character") {
      mplot <- mplot +
        ggtitle(plot_header)
    }

    if (length(IDarray) == 2) {
      mplot <- mplot +
        scale_fill_manual(values = c("#66A53D", "#52307c")) +
        scale_color_manual(values = c("#66A53D", "#52307c"))
      # scale_fill_manual(values = c("#52307c", "#BC6C45")) +
      # scale_color_manual(values = c("#52307c", "#BC6C45"))
      # scale_fill_manual(values = c("blue", "#52307c")) +
      # scale_color_manual(values = c("blue", "#52307c"))
    } else if (length(IDarray) == 3) {
      mplot <- mplot +
        # scale_fill_manual(values = c("#BC6C45", "#66A53D", "#52307c")) +
        # scale_color_manual(values = c("#BC6C45", "#66A53D", "#52307c"))
        scale_fill_manual(values = c("blue", "orange", "#52307c")) +
        scale_color_manual(values = c("blue", "orange", "#52307c"))
    } else {
      mplot <- mplot +
        scale_color_viridis(end = 0.8, discrete = T) +
        scale_fill_viridis(end = 0.8, discrete = T)
    }
    if (spbmsd == T) {
      mplot <- mplot +
        geom_line(
          data = data.frame(xs = spb_float_x, ys = spb_float_y),
          aes(x = xs, y = ys),
          size = 2,
          alpha = 0.5,
          linetype = "dashed"
        )
    }
    if (plot_legend == T) {
      mplot <- mplot +
        theme(legend.position = "bottom")
    } else {
      mplot <- mplot +
        theme(legend.position = "none")
    }




    if (save_plot == T) {
      # Save and Print the plot
      ggsave(
        plot_output_file,
        plot = mplot,
        width = 6,
        height = 4.8
      )
    }

    print(mplot)
    print("Plot's done!")
    setwd(initial_dir)
    global_mmsd <<- msdcomplete
    print(msdaves)
  }
