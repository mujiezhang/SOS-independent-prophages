df1=read.table('prophagedb-qualify_vs-unqualified.tsv',head=TRUE,sep='\t')
library(ggplot2)
library(patchwork) # 用于组合多个ggplot图形
library(ggsignif)

df1$type <- factor(df1$type, levels = c("qualified", "unqualified", "prophagedb"))

mycomparisons2 <- list(c("qualified",'unqualified'),c('unqualified','prophagedb'),c("qualified",'prophagedb'))
p1 <- ggplot(df1, aes(x = type,y=completeness,fill=type)) +
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values=c('qualified'='#F8D186','unqualified'='#E6B081','prophagedb'='#D2D6D7'))+
  geom_signif(
    comparisons = mycomparisons2,
    test = "t.test",  
    step_increase = 0.1,
    map_signif_level = FALSE,
    color = "black",
    size = 0.5
  ) +
  #scale_fill_manual(values=c('#ff9900','#146eb4'))+
  ylab('Completeness (%)')+
  scale_y_continuous(limits = c(0, 130), breaks = seq(0, 130, 20))+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),legend.position ='none')
p1
p2 <- ggplot(df1, aes(x = completeness,fill=type)) +
  geom_density(alpha=0.5) +
  scale_fill_manual(values=c('qualified'='#F8D186','unqualified'='#E6B081','prophagedb'='#D2D6D7'))+
  ylab('Density')+xlab('Completeness (%)')+
  #scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20))+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),legend.position = c(0.6,0.8))
p2
p3 <- ggplot(df1, aes(x = type,y=length/1000,fill=type)) +
  geom_boxplot(outlier.shape = NA)+
  geom_signif(
    comparisons = mycomparisons2,
    test = "t.test",  
    step_increase = 0.1,
    map_signif_level = FALSE,
    color = "black",
    size = 0.5,
    #y_position = 80
  ) +
  scale_fill_manual(values=c('qualified'='#F8D186','unqualified'='#E6B081','prophagedb'='#D2D6D7'))+
  ylab('Genome size (kb)')+ylim(0,100)+
  #scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20))+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),legend.position = 'none')
p3
p4 <- ggplot(df1, aes(x = length/1000,fill=type)) +
  geom_density(alpha=0.5) +
  scale_fill_manual(values=c('qualified'='#F8D186','unqualified'='#E6B081','prophagedb'='#D2D6D7'))+
  ylab('Density')+xlab('Genome size (kb)')+xlim(0,100)+
  #scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20))+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),legend.position = c(0.6,0.8))
p4
layout <- "
ABB
CDD
"

p=p1 + p2 + p3+p4 + plot_layout(design = layout)
p
ggsave(file="length-completeness-quality.pdf",device = cairo_pdf, width=6, height=5,p)

