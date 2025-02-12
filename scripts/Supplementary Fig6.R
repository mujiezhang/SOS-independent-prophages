library(ggplot2)
library(ggsignif)
library(patchwork)

df1=read.table('Supplementary Fig 6-dissimilarity.txt',head=TRUE,sep='\t')
df1$Category <- factor(df1$Category,levels =  c('SiPs','SdPs'))
p1 <- ggplot(df1, aes(x = Category, y = d2star, fill = Category)) + 
  geom_violin(alpha = 1, width = 0.6, trim = TRUE, draw_quantiles = c(0.25, 0.5, 0.75)) +  
  scale_fill_manual(values=c('SiPs'="#ECB477",'SdPs'= "#accde8")) +  
  xlab(NULL)+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   panel.border = element_rect(color = "black", fill = NA, size = 0.5), 
                   axis.text.x = element_text(angle = 0, hjust = 1), 
                   axis.ticks.length.y = unit(0.05, "cm"), 
                   axis.ticks.y = element_line(color = "black"),legend.position = 'none')+
  geom_signif(
    comparisons = list(
      c("SdPs", "SiPs")),
    map_signif_level = F,
    textsize = 3,
    color = "black",
    y_position = max(df1$d2star)*1.05,  
    tip_length = 0,  
    test = "t.test",  
  ) +
  stat_summary(fun = mean, geom = "point", shape = 17, size = 1, color = "black", fill = "black") 

p1
p2 <- ggplot(df1, aes(x = Category, y = d2, fill = Category)) +  
  geom_violin(alpha = 1, width = 0.6, trim = TRUE, draw_quantiles = c(0.25, 0.5, 0.75)) +  
  scale_fill_manual(values=c('SiPs'="#ECB477",'SdPs'= "#accde8")) +  
  xlab(NULL)+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   panel.border = element_rect(color = "black", fill = NA, size = 0.5), 
                   axis.text.x = element_text(angle = 0, hjust = 1), 
                   axis.ticks.length.y = unit(0.05, "cm"), 
                   axis.ticks.y = element_line(color = "black"),legend.position = 'none')+
  geom_signif(
    comparisons = list(
      c("SdPs", "SiPs")),
    map_signif_level = F,
    textsize = 3,
    color = "black",
    y_position = max(df1$d2)*1.05,  
    tip_length = 0,  
    test = "t.test",  
  ) +
  stat_summary(fun = mean, geom = "point", shape = 17, size = 1, color = "black", fill = "black")
p2
p3 <- ggplot(df1, aes(x = Category, y = d2shepp, fill = Category)) +  
  geom_violin(alpha = 1, width = 0.6, trim = TRUE, draw_quantiles = c(0.25, 0.5, 0.75)) +  
  scale_fill_manual(values=c('SiPs'="#ECB477",'SdPs'= "#accde8")) +  
  xlab(NULL)+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   panel.border = element_rect(color = "black", fill = NA, size = 0.5), 
                   axis.text.x = element_text(angle = 0, hjust = 1), 
                   axis.ticks.length.y = unit(0.05, "cm"), 
                   axis.ticks.y = element_line(color = "black"),legend.position = 'none')+
  geom_signif(
    comparisons = list(
      c("SdPs", "SiPs")),
    map_signif_level = F,
    textsize = 3,
    color = "black",
    y_position = max(df1$d2shepp)*1.05,  
    tip_length = 0,  
    test = "t.test",  
  ) +
  stat_summary(fun = mean, geom = "point", shape = 17, size = 1, color = "black", fill = "black")
p3
p4 <- ggplot(df1, aes(x = Category, y = Eu, fill = Category)) +  
  geom_violin(alpha = 1, width = 0.6, trim = TRUE, draw_quantiles = c(0.25, 0.5, 0.75)) + 
  scale_fill_manual(values=c('SiPs'="#ECB477",'SdPs'= "#accde8")) +
  xlab(NULL)+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   panel.border = element_rect(color = "black", fill = NA, size = 0.5), 
                   axis.text.x = element_text(angle = 0, hjust = 1), 
                   axis.ticks.length.y = unit(0.05, "cm"), 
                   axis.ticks.y = element_line(color = "black"),legend.position = 'none')+
  geom_signif(
    comparisons = list(
      c("SdPs", "SiPs")),
    map_signif_level = F,
    textsize = 3,
    color = "black",
    y_position = max(df1$Eu)*1.05, 
    tip_length = 0,  
    test = "t.test",  
  ) +
  stat_summary(fun = mean, geom = "point", shape = 17, size = 1, color = "black", fill = "black")
p4
p5 <- ggplot(df1, aes(x = Category, y = EuF, fill = Category)) + 
  geom_violin(alpha = 1, width = 0.6, trim = TRUE, draw_quantiles = c(0.25, 0.5, 0.75)) +  
  scale_fill_manual(values=c('SiPs'="#ECB477",'SdPs'= "#accde8")) + 
  xlab(NULL)+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   panel.border = element_rect(color = "black", fill = NA, size = 0.5), 
                   axis.text.x = element_text(angle = 0, hjust = 1), 
                   axis.ticks.length.y = unit(0.05, "cm"), 
                   axis.ticks.y = element_line(color = "black"),legend.position = 'none')+
  geom_signif(
    comparisons = list(
      c("SdPs", "SiPs")),
    map_signif_level = F,
    textsize = 3,
    color = "black",
    y_position = max(df1$EuF)*1.05,
    tip_length = 0,  
    test = "t.test", 
  ) +
  stat_summary(fun = mean, geom = "point", shape = 17, size = 1, color = "black", fill = "black")  # 三角形表示平均值
p5
p6 <- ggplot(df1, aes(x = Category, y = Hao, fill = Category)) + 
  geom_violin(alpha = 1, width = 0.6, trim = TRUE, draw_quantiles = c(0.25, 0.5, 0.75)) + 
  scale_fill_manual(values=c('SiPs'="#ECB477",'SdPs'= "#accde8")) +
  xlab(NULL)+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   panel.border = element_rect(color = "black", fill = NA, size = 0.5), 
                   axis.text.x = element_text(angle = 0, hjust = 1), 
                   axis.ticks.length.y = unit(0.05, "cm"), 
                   axis.ticks.y = element_line(color = "black"),legend.position = 'none')+
  geom_signif(
    comparisons = list(
      c("SdPs", "SiPs")),
    map_signif_level = F,
    textsize = 3,
    color = "black",
    y_position = max(df1$Hao)*1.05, 
    tip_length = 0,  
    test = "t.test",  
  ) +
  stat_summary(fun = mean, geom = "point", shape = 17, size = 1, color = "black", fill = "black") 
p6
p7 <- ggplot(df1, aes(x = Category, y = JS, fill = Category)) +  
  geom_violin(alpha = 1, width = 0.6, trim = TRUE, draw_quantiles = c(0.25, 0.5, 0.75)) +  
  scale_fill_manual(values=c('SiPs'="#ECB477",'SdPs'= "#accde8")) +  
  xlab(NULL)+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   panel.border = element_rect(color = "black", fill = NA, size = 0.5), 
                   axis.text.x = element_text(angle = 0, hjust = 1), 
                   axis.ticks.length.y = unit(0.05, "cm"), 
                   axis.ticks.y = element_line(color = "black"),legend.position = 'none')+
  geom_signif(
    comparisons = list(
      c("SdPs", "SiPs")),
    map_signif_level = F,
    textsize = 3,
    color = "black",
    y_position = max(df1$JS)*1.05,  
    tip_length = 0,  
    test = "t.test",  
  ) +
  stat_summary(fun = mean, geom = "point", shape = 17, size = 1, color = "black", fill = "black") 
p7
p8 <- ggplot(df1, aes(x = Category, y = Ma, fill = Category)) + 
  geom_violin(alpha = 1, width = 0.6, trim = TRUE, draw_quantiles = c(0.25, 0.5, 0.75)) + 
  scale_fill_manual(values=c('SiPs'="#ECB477",'SdPs'= "#accde8")) + 
  xlab(NULL)+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   panel.border = element_rect(color = "black", fill = NA, size = 0.5), 
                   axis.text.x = element_text(angle = 0, hjust = 1), 
                   axis.ticks.length.y = unit(0.05, "cm"), 
                   axis.ticks.y = element_line(color = "black"),legend.position = 'none')+
  geom_signif(
    comparisons = list(
      c("SdPs", "SiPs")),
    map_signif_level = F,
    textsize = 3,
    color = "black",
    y_position = max(df1$Ma)*1.05,
    tip_length = 0,  
    test = "t.test", 
  ) +
  stat_summary(fun = mean, geom = "point", shape = 17, size = 1, color = "black", fill = "black") 
p8
p9 <- ggplot(df1, aes(x = Category, y = Ch, fill = Category)) + 
  geom_violin(alpha = 1, width = 0.6, trim = TRUE, draw_quantiles = c(0.25, 0.5, 0.75)) + 
  scale_fill_manual(values=c('SiPs'="#ECB477",'SdPs'= "#accde8")) + 
  xlab(NULL)+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   panel.border = element_rect(color = "black", fill = NA, size = 0.5), 
                   axis.text.x = element_text(angle = 0, hjust = 1), 
                   axis.ticks.length.y = unit(0.05, "cm"), 
                   axis.ticks.y = element_line(color = "black"),legend.position = 'none')+
  geom_signif(
    comparisons = list(
      c("SdPs", "SiPs")),
    map_signif_level = F,
    textsize = 3,
    color = "black",
    y_position = max(df1$Ch)*1.05,
    tip_length = 0, 
    test = "t.test", 
  ) +
  stat_summary(fun = mean, geom = "point", shape = 17, size = 1, color = "black", fill = "black")
p9
design1 <- "
  1234
  5678
  "
ppp=p2+p3+p4+p5+p6+p7+p8+p9+ plot_layout(design = design1)
ppp

