# code adapted from https://forum.qiime2.org/t/tutorial-integrating-qiime2-and-r-for-data-visualization-and-analysis-using-qiime2r/4121
install.packages("tidyverse")
library(tidyverse)

depths <- c(1, 5, 10, 30, 60)

for (depth in depths) {
  # Generate file paths
  input_file <- paste0('R_testing_vis/', depth, '_chloroplastfor_R.csv')
  outputpath <- paste0('R_plots/PCoA_chloroplast_', depth, 'm.png')
  
  # Read the data
  df <- read.csv(input_file)
  
  df$size_code <- factor(df$size_code, 
                         levels = c("L", "S", "SL", "W"),
                         labels = c("Large", "Small", "Small+Large", "Whole"))
  

  gg<- ggplot(df, aes(x=dim1, y=dim2, color=size_code, shape=Time))+
    scale_color_manual(values = c("#d55e00", "#5975a4","#b55d60", "#5f9e6e"), name ='Size fraction') +#, size=weekn)) +
    geom_point(alpha=0.5, size=4)+ #alpha controls transparency and helps when points are overlapping
    scale_shape_manual(values=c(16,1), name="Time") + #see http://www.sthda.com/sthda/RDoc/figure/graphs/r-plot-pch-symbols-points-in-r.png for numeric shape codes
    #scale_size_continuous(name="Shannon Diversity") +
    #scale_color_discrete(name="Size fraction") +
    theme_minimal() +
    labs(
      x = paste0("PCo 1 (", df$PCo.1[1], "%)"),
      y = paste0("PCo 2 (", df$PCo.2[1], "%)")
    ) +
    stat_ellipse(data = df, aes(x=dim1, y=dim2,color=size_code, group=size_code),type = "norm", size = 1.2)+
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=18,face="bold"))
  plot(gg)
  
  gg_no_legend <- gg + theme(legend.position = "none")
  ggsave(filename=outputpath, gg, width = 6, height = 6, dpi = 300, scale = 0.5)
}
#ggsave(filename = outputpath, width = 5, height = 4, device='tiff', dpi=700) # save a PDF 3 inches by 4 inches

