multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

plot_error_and_variance <- function(input_file){
  library('ggplot2')
  load(input_file)
  total_nsig = (nsig_max - nsig_min + 1)*(nsig_min + nsig_max)/2
  print(total_nsig)  

  signature_number_column = rep(0, total_nsig)
  error_column = rep(0, total_nsig)

  count = 1
  for(isig in nsig_min:nsig_max){ 
    for(i in 1:isig){
      signature_number_column[count] = isig
      error_column[count] = min_error[isig]
      count = count + 1
    }
  }
  print(signature_number_column)
#  png(sprintf('bars_num%d_totalnum%d.png',k,isig), width=6000, height=1500, res=300)
  frame_error_variance <- data.frame(error_column, variance_all, signature_number_column)
  colnames(frame_error_variance) = c('error', 'variance', 'signature')
  variance_plot <- ggplot(frame_error_variance, aes(x = signature, y = variance)) +  geom_boxplot(aes(group = cut_width(signature, 1)),outlier.size=2,outlier.colour="green") 
  error_plot <- ggplot(frame_error_variance, aes(x = signature, y = error)) +  geom_boxplot(aes(group = cut_width(signature, 1)),outlier.size=2,outlier.colour="green") 
#+  stat_summary(fun.y="mean", geom = "point", shape=23, size =3, fill="pink") + 
#  ggtitle("Portugese Sea Battles")  
  multiplot(variance_plot, error_plot, cols = 1)
  ggsave(filename = sprintf("plots_error_variance/%s.png",input_file))

  
}