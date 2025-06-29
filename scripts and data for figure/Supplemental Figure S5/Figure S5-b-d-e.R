library(ggplot2)
library(ggsignif)
library(ggThemeAssist)
library(patchwork)

#for bar plot
df<- read.table("class-distribution.tsv",header = T,sep='\t')
df$class<- factor(df$class, levels = c('Gammaproteobacteria', 'Bacilli', 'Actinomycetia', 'Clostridia', 'Alphaproteobacteria', 
                                       'Spirochaetia', 'Bacteroidia', 'Planctomycetia', 
                                       'Thermoanaerobacteria', 'Negativicutes', 'Desulfobacteria', 'Dehalobacteriia', 
                                       'Cyanobacteriia', 'Coriobacteriia','Unclassified Bacteria'))

p=ggplot(df, aes( x = class,y=gemomes))+
  geom_col(width = 0.7,colour='black',linewidth=0.2,fill='#adcde8')+
  ylim(0,18000)+ylab('Genomes')+xlab(NULL)+
  geom_text(aes(label = gemomes), vjust = -1.5)+
  theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   legend.position = 'none')
p
ggsave(file="barplot.pdf",device = cairo_pdf, width=6, height=4,p)

##for violin
df1=read.table("2-cluster.txt",head=TRUE,sep='\t')
df2=read.table("3-cluster.txt",head=TRUE,sep='\t')
df3=read.table("1-cluster.txt",head=TRUE,sep='\t')
df4=read.table("4-cluster.txt",head=TRUE,sep='\t')
df5=read.table("5-cluster.txt",head=TRUE,sep='\t')
df6=read.table("6-cluster.txt",head=TRUE,sep='\t')

df1$Labels <- as.factor(df1$Labels)
df2$Labels <- as.factor(df2$Labels)
df3$Labels <- as.factor(df3$Labels)
df4$Labels <- as.factor(df4$Labels)
df5$Labels <- as.factor(df5$Labels)
df6$Labels <- as.factor(df6$Labels)
p1 <- ggplot() +
  geom_violin(df1, mapping=aes(x=name, y=HI),trim=FALSE, fill='#d6dce5', alpha=0) +
  xlab(NULL) +
  ylab('Heterology Index') +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = 'none') +
  geom_jitter(df1, mapping=aes(x=name, y=HI,fill=Labels),width = 0.08,size = 1.5,alpha = 0.5,shape=21)
p1

p2 <- ggplot() +
  geom_violin(df2, mapping=aes(x=name, y=HI),trim=FALSE, fill='#d6dce5', alpha=0) +
  xlab(NULL) +
  ylab('Heterology Index') +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = 'none') +
  geom_jitter(df2, mapping=aes(x=name, y=HI,fill=Labels),width = 0.08,size = 1.5,alpha = 0.5,shape=21)
p2
p3 <- ggplot() +
  geom_violin(df3, mapping=aes(x=name, y=HI),trim=FALSE, fill='#d6dce5', alpha=0) +
  xlab(NULL) +
  ylab('Heterology Index') +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = 'none') +
  geom_jitter(df3, mapping=aes(x=name, y=HI,fill=Labels),width = 0.08,size = 1.5,alpha = 0.5,shape=21)
p3
p4 <- ggplot() +
  geom_violin(df4, mapping=aes(x=name, y=HI),trim=FALSE, fill='#d6dce5', alpha=0) +
  xlab(NULL) +
  ylab('Heterology Index') +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = 'none') +
  geom_jitter(df4, mapping=aes(x=name, y=HI,fill=Labels),width = 0.08,size = 1.5,alpha = 0.5,shape=21)
p4
p5 <- ggplot() +
  geom_violin(df5, mapping=aes(x=name, y=HI),trim=FALSE, fill='#d6dce5', alpha=0) +
  xlab(NULL) +
  ylab('Heterology Index') +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = 'none') +
  geom_jitter(df5, mapping=aes(x=name, y=HI,fill=Labels),width = 0.08,size = 1.5,alpha = 0.5,shape=21)
p5
p6 <- ggplot() +
  geom_violin(df6, mapping=aes(x=name, y=HI),trim=FALSE, fill='#d6dce5', alpha=0) +
  xlab(NULL) +
  ylab('Heterology Index') +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = 'none') +
  geom_jitter(df6, mapping=aes(x=name, y=HI,fill=Labels),width = 0.08,size = 1.5,alpha = 0.5,shape=21)
p6

pp=p3+p1+p2+p4+p5+p6
pp
ggsave(file="1-2-3-4-5-6-clusters.pdf",device = cairo_pdf, width=9, height=7,pp)



#for boxplot
df7<- read.table("3-type-hiscore.tsv",header = T)
#df$type <- factor(df$type, levels = c("re", "po", "no","repo", "reno", "pono"))
df1$type<- factor(df1$type, levels = c('sip','sup','sdp'))
mycomparisons1 <- list(c("sdp",'sup'),c("sdp", "sip"),c("sup",'sip'))
p7=ggplot(df7, aes( x = type,y=hi,fill=type))+
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values=c('re'='#accde8','po'='#ecf4dd','no'='#ecb477'))+
  labs(x = NULL,
       y = "score") +
  geom_signif(
    comparisons = mycomparisons1,
    test =t.test,
    y_position = c(22,25,28),
    tip_length = 0.01,
  )+
  geom_jitter(width=0.1,size=0.5,fill='black')+
  theme_bw()+theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = 'none')
p7
ggsave(file="boxplot.pdf",device = cairo_pdf, width=5, height=5,p7)




