df1=read.table('ecoli-k12-HI-predict-all.tsv',head=TRUE,sep='\t')
df2=read.table('25point-all.tsv',head=TRUE,sep='\t')

library(ggThemeAssist)
library(ggplot2) 
library(patchwork)

a<-subset(df1, model == "MeanShift")
aa<-subset(df2, model == "MeanShift")
p1 <- ggplot() +
  geom_violin(a, mapping=aes(x=model, y=hi),trim=FALSE, fill='#d6dce5', alpha=0) +
  geom_jitter(a, mapping=aes(x=model, y=hi,fill='grey',color='grey'),width = 0.08,size = 1,alpha = 0.5,shape=21)+
  xlab(NULL) +
  ylab('Heterology Index') +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = 'none') +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 13.20, fill="#accde8", alpha=0.3) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 13.20, ymax = 14.78, fill="#ECF4DD",  alpha=0.3) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 14.78, ymax = Inf, fill="#ECB477", alpha=0.3) +
  geom_jitter(aa, mapping=aes(x=model, y=hi,fill=label),width = 0.08,size = 1.5,alpha = 1,shape=22)+
  scale_fill_manual(values = c('yes'="#accde8",'no'="#ECB477")) 
p1

b<-subset(df1, model == "KMeans")
bb<-subset(df2, model == "KMeans")
p2 <- ggplot() +
  geom_violin(b, mapping=aes(x=model, y=hi),trim=FALSE, fill='#d6dce5', alpha=0) +
  geom_jitter(b, mapping=aes(x=model, y=hi,fill='grey',color='grey'),width = 0.08,size = 1,alpha = 0.5,shape=21)+
  xlab(NULL) +
  ylab('Heterology Index') +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = 'none') +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 17.9124, fill="#accde8", alpha=0.3) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 17.9124, ymax = 19.48, fill="#ECF4DD",  alpha=0.3) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 19.48, ymax = Inf, fill="#ECB477", alpha=0.3) +
  geom_jitter(bb, mapping=aes(x=model, y=hi,fill=label),width = 0.08,size = 1.5,alpha = 1,shape=22)+
  scale_fill_manual(values = c('yes'="#accde8",'no'="#ECB477")) 
p2
c<-subset(df1, model == "Birch")
cc<-subset(df2, model == "Birch")
p3 <- ggplot() +
  geom_violin(c, mapping=aes(x=model, y=hi),trim=FALSE, fill='#d6dce5', alpha=0) +
  geom_jitter(c, mapping=aes(x=model, y=hi,fill='grey',color='grey'),width = 0.08,size = 1,alpha = 0.5,shape=21)+
  xlab(NULL) +
  ylab('Heterology Index') +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = 'none') +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 19.75898, fill="#accde8", alpha=0.3) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 19.75898, ymax = 21.330265, fill="#ECF4DD",  alpha=0.3) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin =21.330265, ymax = Inf, fill="#ECB477", alpha=0.3) +
  geom_jitter(cc, mapping=aes(x=model, y=hi,fill=label),width = 0.08,size = 1.5,alpha = 1,shape=22)+
  scale_fill_manual(values = c('yes'="#accde8",'no'="#ECB477")) 
p3
d<-subset(df1, model == "GaussianMixture")
dd<-subset(df2, model == "GaussianMixture")
p4 <- ggplot() +
  geom_violin(d, mapping=aes(x=model, y=hi),trim=FALSE, fill='#d6dce5', alpha=0) +
  geom_jitter(d, mapping=aes(x=model, y=hi,fill='grey',color='grey'),width = 0.08,size = 1,alpha = 0.5,shape=21)+
  xlab(NULL) +
  ylab('Heterology Index') +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = 'none') +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 15.7386275, fill="#accde8", alpha=0.3) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 15.7386275, ymax = 17.30991, fill="#ECF4DD",  alpha=0.3) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin =17.30991, ymax = Inf, fill="#ECB477", alpha=0.3) +
  geom_jitter(dd, mapping=aes(x=model, y=hi,fill=label),width = 0.08,size = 1.5,alpha = 1,shape=22)+
  scale_fill_manual(values = c('yes'="#accde8",'no'="#ECB477")) 
p4
e<-subset(df1, model == "BayesianGaussianMixture")
ee<-subset(df2, model == "BayesianGaussianMixture")
p5 <- ggplot() +
  geom_violin(e, mapping=aes(x=model, y=hi),trim=FALSE, fill='#d6dce5', alpha=0) +
  geom_jitter(e, mapping=aes(x=model, y=hi,fill='grey',color='grey'),width = 0.08,size = 1,alpha = 0.5,shape=21)+
  xlab(NULL) +
  ylab('Heterology Index') +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = 'none') +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 11.87721, fill="#accde8", alpha=0.3) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 11.87721, ymax = 13.448495, fill="#ECF4DD",  alpha=0.3) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin =13.448495, ymax = Inf, fill="#ECB477", alpha=0.3) +
  geom_jitter(ee, mapping=aes(x=model, y=hi,fill=label),width = 0.08,size = 1.5,alpha = 1,shape=22)+
  scale_fill_manual(values = c('yes'="#accde8",'no'="#ECB477")) 
p5
f<-subset(df1, model == "AgglomerativeClustering")
ff<-subset(df2, model == "AgglomerativeClustering")
p6 <- ggplot() +
  geom_violin(f, mapping=aes(x=model, y=hi),trim=FALSE, fill='#d6dce5', alpha=0) +
  geom_jitter(f, mapping=aes(x=model, y=hi,fill='grey',color='grey'),width = 0.08,size = 1,alpha = 0.5,shape=21)+
  xlab(NULL) +
  ylab('Heterology Index') +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = 'none') +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 20.1504275, fill="#accde8", alpha=0.3) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 20.1504275, ymax = 21.721712500000002, fill="#ECF4DD",  alpha=0.3) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin =21.721712500000002, ymax = Inf, fill="#ECB477", alpha=0.3) +
  geom_jitter(ff, mapping=aes(x=model, y=hi,fill=label),width = 0.08,size = 1.5,alpha = 1,shape=22)+
  scale_fill_manual(values = c('yes'="#accde8",'no'="#ECB477")) 
p6
g<-subset(df1, model == "SpectralClustering")
gg<-subset(df2, model == "SpectralClustering")
p7 <- ggplot() +
  geom_violin(g, mapping=aes(x=model, y=hi),trim=FALSE, fill='#d6dce5', alpha=0) +
  geom_jitter(g, mapping=aes(x=model, y=hi,fill='grey',color='grey'),width = 0.08,size = 1,alpha = 0.5,shape=21)+
  xlab(NULL) +
  ylab('Heterology Index') +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = 'none') +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax =10.808887499999999, fill="#accde8", alpha=0.3) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 10.808887499999999, ymax = 12.3801725, fill="#ECF4DD",  alpha=0.3) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin =12.3801725, ymax = Inf, fill="#ECB477", alpha=0.3) +
  geom_jitter(gg, mapping=aes(x=model, y=hi,fill=label),width = 0.08,size = 1.5,alpha = 1,shape=22)+
  scale_fill_manual(values = c('yes'="#accde8",'no'="#ECB477")) 
p7
layout <- "
ABCD
EFG#
"

p=p1+p2+p3+p4+p5+p6+p7 + plot_layout(design = layout)
p
ggsave(file="model-compare.pdf",device = cairo_pdf, width=10, height=5,p)