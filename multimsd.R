



mmsdplot <-
  function(IDarray,
           spbmsd = T,
           title = "test",
           plot_header = NA,
           merge_input_dir = paste(dataset_dir, "trajectory_r_data/analytical_environment/MSD_merged", sep = ""),
           plot_output_dir,
           plot_legend = F,
           save_plot = F) {
    setwd(merge_input_dir)
    msdcomplete <-
      data.frame(
        MSD = numeric(),
        taus = numeric(),
        experiment = character(),
        ID = character()
      )
    for (i in 1:length(IDarray)) {
      fname <- list.files(pattern = IDarray[i])
      print(fname)
      dffloat <- fread(list.files(pattern = IDarray[i]))
      # print(dffloat[,c("MSD","taus","experiment")])
      msdcomplete <-
        rbind(msdcomplete, cbind(
          dffloat[, c("MSD", "taus", "experiment")],
          data.frame(ID = as.character(IDarray[i]))
        ))
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
      geom_smooth(data = msdcomplete, aes(
        x = taus,
        y = MSD,
        color = ID,
        fill = ID
      )) +
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
        theme(legend.position = "right")
    } else {
      mplot <- mplot +
        theme(legend.position = "none")
    }




    if (save_plot == T) {
      # Save and Print the plot
      ggsave(
        plot_output_dir,
        plot = mplot,
        width = 6,
        height = 4.8
      )
    }

    print(mplot)
    print("Plot's done!")

    global_mmsd <<- msdcomplete
  }
