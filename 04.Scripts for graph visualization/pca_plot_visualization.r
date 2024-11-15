library(ggplot2)
library(gridExtra)
mydata1 <-  read.table("NEW111",header = TRUE)
m_colour1=c("Ancient"="#a4a4a4","ED"="#9A50D9","EWB"="#FF0000",
           "SC"="#2AB359","NCWB"="#0000C8","SCWB"="#FFAD00","NC"="#00C8F0")
p1 <- ggplot(mydata1, aes(x=PC1,y=PC2, color=group)) + geom_point(size=3) + 
  scale_color_manual(values=m_colour1,breaks=c("ED","EWB","NCWB","NC","SCWB","SC","Ancient")) +
  theme_bw() + 
  xlab("PC1(31.15%)") + 
  ylab("PC2(15.37%)") + 
  guides(colour = guide_legend(override.aes = list(shape = 16))) +
  theme(panel.grid = element_blank()) +
  labs(title = "Auto")
p1
