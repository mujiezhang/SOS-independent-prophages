
sj=read.table('Serratia_J.faa_intergenic_region.fna_pattern.txt_whole_genome_HI.tsv',head=TRUE,sep='\t')
sj2=read.table('Serratia_J_point.tsv',head=TRUE,sep='\t')
ha=read.table('Hafnia.faa_intergenic_region.fna_pattern.txt_whole_genome_HI.tsv',head=TRUE,sep='\t')
ha2=read.table('Hafnia_point.tsv',head=TRUE,sep='\t')
LT2=read.table('Salmonella-LT2-Fels1-host.faa_intergenic_region.fna_pattern.txt_whole_genome_HI.txt',head=TRUE,sep='\t')
LT22=read.table('Salmonella-LT2.tsv',head=TRUE,sep='\t')
vibrio=read.table('Vibrio_K01M1.faa_intergenic_region.fna_pattern.txt_whole_genome_HI.txt',head=TRUE,sep='\t')
vibrio2=read.table('Vibrio_K01M1.tsv',head=TRUE,sep='\t')
hs=read.table('ecoli-HS.faa_intergenic_region.fna_pattern.txt_whole_genome_HI.txt',head=TRUE,sep='\t')
hs2=read.table('ecoli-HS_point.tsv',head=TRUE,sep='\t')
st8624=read.table('ecoli-8624.faa_intergenic_region.fna_pattern.txt_whole_genome_HI.txt',head=TRUE,sep='\t')
st86242=read.table('ecoli-8624_point.tsv',head=TRUE,sep='\t')
pao1=read.table('pao1.faa_intergenic_region.fna_pattern.txt_whole_genome_HI.txt',head=TRUE,sep='\t')
pao12=read.table('pao1_point.tsv',head=TRUE,sep='\t')
library(ggThemeAssist)
library(ggplot2) 
library(patchwork)  

p1<- ggplot(sj,aes(x=name,y=HI))+
  geom_violin(trim=FALSE,fill='#d6dce5',alpha=0.5)+xlab(NULL)+ylab('Heterology Index')+
  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = 'none')+
  
  geom_jitter(sj,mapping=aes(x=name,y=HI),width=0.1,size=3,alpha=0.5, color='grey',fill='grey')+
  geom_jitter(sj2,mapping=aes(x=name,y=HI,fill=state),shape=22,width=0.1,size=3,alpha=0.7,color='black')+
  scale_fill_manual(values = c('no'= "#ECB477",'yes'="#ACCDE8"))+
  geom_hline(aes(yintercept=12.71), colour="#FA7F6F", linetype="dashed") +
  geom_hline(aes(yintercept=14.11), colour="#FA7F6F", linetype="dashed") +
  annotate("text", x=0.62, y=15, label="14.11", color="black")+
  annotate("text", x=0.62, y=12, label="12.71", color="black")
p1

p2<- ggplot(ha,aes(x=name,y=HI))+
  geom_violin(trim=FALSE,fill='#d6dce5',alpha=0.5)+xlab(NULL)+ylab('Heterology Index')+
  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = 'none')+
  
  geom_jitter(ha,mapping=aes(x=name,y=HI),width=0.1,size=3,alpha=0.5, color='grey',fill='grey')+
  geom_jitter(ha2,mapping=aes(x=name,y=HI,fill=state),shape=22,width=0.1,size=3,alpha=0.7,color='black')+
  scale_fill_manual(values = c('no'= "#ECB477",'yes'="#ACCDE8"))+
  geom_hline(aes(yintercept=11.85), colour="#FA7F6F", linetype="dashed") +
  geom_hline(aes(yintercept=13.38), colour="#FA7F6F", linetype="dashed") +
  annotate("text", x=0.62, y=13.9, label="13.38", color="black")+
  annotate("text", x=0.62, y=11.3, label="11.85", color="black")
p2

p3<- ggplot(LT2,aes(x=name,y=HI))+
  geom_violin(trim=FALSE,fill='#d6dce5',alpha=0.5)+xlab(NULL)+ylab('Heterology Index')+
  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = 'none')+
  
  geom_jitter(LT2,mapping=aes(x=name,y=HI),width=0.1,size=3,alpha=0.5, color='grey',fill='grey')+
  geom_jitter(LT22,mapping=aes(x=name,y=HI,fill=state),shape=22,width=0.1,size=3,alpha=0.7,color='black')+
  scale_fill_manual(values = c('no'= "#ECB477",'yes'="#ACCDE8"))+
  geom_hline(aes(yintercept=14.99), colour="#FA7F6F", linetype="dashed") +
  geom_hline(aes(yintercept=16.68), colour="#FA7F6F", linetype="dashed") +
  annotate("text", x=0.62, y=14.5, label="14.99", color="black")+
  annotate("text", x=0.62, y=17.2, label="16.68", color="black")
p3
p4<- ggplot(vibrio,aes(x=name,y=HI))+
  geom_violin(trim=FALSE,fill='#d6dce5',alpha=0.5)+xlab(NULL)+ylab('Heterology Index')+
  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = 'none')+
  
  geom_jitter(vibrio,mapping=aes(x=name,y=HI),width=0.1,size=3,alpha=0.5, color='grey',fill='grey')+
  geom_jitter(vibrio2,mapping=aes(x=name,y=HI,fill=state),shape=22,width=0.1,size=3,alpha=0.7,color='black')+
  scale_fill_manual(values = c('no'= "#ECB477",'yes'="#ACCDE8"))+
  geom_hline(aes(yintercept=15.49), colour="#FA7F6F", linetype="dashed") +
  geom_hline(aes(yintercept=16.84), colour="#FA7F6F", linetype="dashed") +
  annotate("text", x=0.62, y=15, label="15.49", color="black")+
  annotate("text", x=0.62, y=17.3, label="16.84", color="black")
p4
p5<- ggplot(hs,aes(x=name,y=HI))+
  geom_violin(trim=FALSE,fill='#d6dce5',alpha=0.5)+xlab(NULL)+ylab('Heterology Index')+
  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = 'none')+
  
  geom_jitter(hs,mapping=aes(x=name,y=HI),width=0.1,size=3,alpha=0.5, color='grey',fill='grey')+
  geom_jitter(hs2,mapping=aes(x=name,y=HI,fill=state),shape=22,width=0.1,size=3,alpha=0.7,color='black')+
  scale_fill_manual(values = c('no'= "#ECB477",'yes'="#ACCDE8"))+
  geom_hline(aes(yintercept=14.28), colour="#FA7F6F", linetype="dashed") +
  geom_hline(aes(yintercept=15.85), colour="#FA7F6F", linetype="dashed") +
  annotate("text", x=0.62, y=13.70, label="14.28", color="black")+
  annotate("text", x=0.62, y=16.35, label="15.85", color="black")
p5
p6<- ggplot(st8624,aes(x=name,y=HI))+
  geom_violin(trim=FALSE,fill='#d6dce5',alpha=0.5)+xlab(NULL)+ylab('Heterology Index')+
  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = 'none')+
  
  geom_jitter(st8624,mapping=aes(x=name,y=HI),width=0.1,size=3,alpha=0.5, color='grey',fill='grey')+
  geom_jitter(st86242,mapping=aes(x=name,y=HI,fill=state),shape=22,width=0.1,size=3,alpha=0.7,color='black')+
  scale_fill_manual(values = c('no'= "#ECB477",'yes'="#ACCDE8"))+
  geom_hline(aes(yintercept=15.06), colour="#FA7F6F", linetype="dashed") +
  geom_hline(aes(yintercept=16.61), colour="#FA7F6F", linetype="dashed") +
  annotate("text", x=0.62, y=14.5, label="15.06", color="black")+
  annotate("text", x=0.62, y=17.1, label="16.61", color="black")
p6
p7<- ggplot(pao1,aes(x=name,y=HI))+
  geom_violin(trim=FALSE,fill='#d6dce5',alpha=0.5)+xlab(NULL)+ylab('Heterology Index')+
  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = 'none')+
  
  geom_jitter(pao1,mapping=aes(x=name,y=HI),width=0.1,size=3,alpha=0.5, color='grey',fill='grey')+
  geom_jitter(pao12,mapping=aes(x=name,y=HI,fill=state),shape=22,width=0.1,size=3,alpha=0.7,color='black')+
  scale_fill_manual(values = c('no'= "#ECB477",'yes'="#ACCDE8"))+
  geom_hline(aes(yintercept=16.88), colour="#FA7F6F", linetype="dashed") +
  geom_hline(aes(yintercept=18.59), colour="#FA7F6F", linetype="dashed") +
  annotate("text", x=0.62, y=16.38, label="16.88", color="black")+
  annotate("text", x=0.62, y=19, label="18.59", color="black")
p7
p=p1/p2/p3/p4/p5/p6/p7
p
ggsave(file="all_positive.pdf",device = cairo_pdf, width=4, height=24,p)
