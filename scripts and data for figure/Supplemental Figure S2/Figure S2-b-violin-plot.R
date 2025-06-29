library(ggThemeAssist)
library(ggplot2) 
library(patchwork)  
wp3=read.table('wp3.faa_intergenic_region.fna_pattern.txt_whole_genome_HI.txt',head=TRUE,sep='\t')
wp33=read.table('wp3_point.tsv',head=TRUE,sep='\t')
p1<- ggplot(wp3,aes(x=name,y=HI))+
  geom_violin(trim=FALSE,fill='#d6dce5',alpha=0.5)+xlab(NULL)+ylab('Heterology Index')+
  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = 'none')+
  geom_jitter(wp3,mapping=aes(x=name,y=HI),width=0.1,size=3,alpha=0.5, color='grey',fill='grey')+
  geom_jitter(wp33,mapping=aes(x=name,y=HI,fill=state),shape=22,width=0.1,size=3,alpha=0.7,color='black')+
  scale_fill_manual(values = c('no'= "#ECB477",'yes'="#ACCDE8"))+
  geom_hline(aes(yintercept=12.71), colour="#FA7F6F", linetype="dashed") +
  geom_hline(aes(yintercept=14.11), colour="#FA7F6F", linetype="dashed") +
  annotate("text", x=0.62, y=15, label="14.11", color="black")+
  annotate("text", x=0.62, y=12, label="12.71", color="black")
p1
