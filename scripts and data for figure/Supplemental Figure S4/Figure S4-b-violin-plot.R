library(ggThemeAssist)
library(ggplot2) 
library(patchwork) 

wp2=read.table('wp2.faa_intergenic_region.fna_pattern.txt_whole_genome_HI.txt',head=TRUE,sep='\t')
wp22=read.table('wp2_point.tsv',head=TRUE,sep='\t')

p1<- ggplot(wp2,aes(x=name,y=HI))+
  geom_violin(trim=FALSE,fill='#d6dce5',alpha=0.5)+xlab(NULL)+ylab('Heterology Index')+
  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = 'none')+
  
  geom_jitter(wp2,mapping=aes(x=name,y=HI),width=0.1,size=3,alpha=0.5, color='grey',fill='grey')+
  geom_jitter(wp22,mapping=aes(x=name,y=HI,fill=state),shape=22,width=0.1,size=3,alpha=0.7,color='black')+
  scale_fill_manual(values = c('no'= "#ECB477",'yes'="#ACCDE8"))+
  geom_hline(aes(yintercept=11.85), colour="#FA7F6F", linetype="dashed") +
  geom_hline(aes(yintercept=13.38), colour="#FA7F6F", linetype="dashed") +
  annotate("text", x=0.62, y=13.9, label="13.38", color="black")+
  annotate("text", x=0.62, y=11.3, label="11.85", color="black")
p1