# code adapted from https://forum.qiime2.org/t/tutorial-integrating-qiime2-and-r-for-data-visualization-and-analysis-using-qiime2r/4121
library(tidyverse)

df<-read_csv('R_testing_vis/1chloroplastfor_R.csv')
outputpath<-"R_testing_vis/plots/PCoA_chloro_1m.png"
head(df)

gg<- ggplot(df, aes(x=dim1, y=dim2, color=size_code, shape=Time))+
  scale_color_manual(values = c("#d55e00", "#5975a4","#b55d60", "#5f9e6e"), name ='Size fraction') +#, size=weekn)) +
  geom_point(alpha=0.5) + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  scale_shape_manual(values=c(16,1), name="Time") + #see http://www.sthda.com/sthda/RDoc/figure/graphs/r-plot-pch-symbols-points-in-r.png for numeric shape codes
  #scale_size_continuous(name="Shannon Diversity") +
  #scale_color_discrete(name="Size fraction") +
  labs(y="PCo 2 (27%)", 
       x="PCo 1 (42%)") +
  stat_ellipse(data = df, aes(x=dim1, y=dim2,color=size_code, group=size_code),type = "norm")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))
plot(gg)
ggsave(filename = outputpath, width = 5, height = 4, device='tiff', dpi=700) # save a PDF 3 inches by 4 inches

