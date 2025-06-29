df1=read.table('host-HI-infor-top10genus-with-family.tsv',head=TRUE,sep='\t')
df2=read.table('host-HI-infor-top10genus-with-family-ave.tsv',head=TRUE,sep='\t')
df11=read.table('host-HI-infor-top10species-with-genus.tsv',head=TRUE,sep='\t')
df22=read.table('host-HI-infor-top10species-with-genus-ave.tsv',head=TRUE,sep='\t')
library(ggplot2)  
library(patchwork)  


df2$genus<- factor(df2$genus,levels = c('g__Escherichia','g__Klebsiella','g__Salmonella','g__Enterobacter','g__Citrobacter','g__Pseudomonas',
                                        'g__Pseudomonas_E','g__Vibrio','g__Burkholderia','g__Aeromonas'))
df1$genus<- factor(df1$genus,levels =c('g__Escherichia','g__Klebsiella','g__Salmonella','g__Enterobacter','g__Citrobacter','g__Pseudomonas',
                                       'g__Pseudomonas_E','g__Vibrio','g__Burkholderia','g__Aeromonas'))
df22$species<- factor(df22$species,levels = c('s__Escherichia coli', 's__Klebsiella pneumoniae',
                                             's__Klebsiella variicola','s__Klebsiella quasipneumoniae','s__Salmonella enterica',
                                             's__Pseudomonas aeruginosa','s__Enterobacter hormaechei_A',
                                             's__Citrobacter freundii','s__Proteus mirabilis','s__Burkholderia mallei'))
df11$species<- factor(df11$species,levels = c('s__Escherichia coli', 's__Klebsiella pneumoniae',
                                             's__Klebsiella variicola','s__Klebsiella quasipneumoniae','s__Salmonella enterica',
                                             's__Pseudomonas aeruginosa','s__Enterobacter hormaechei_A',
                                             's__Citrobacter freundii','s__Proteus mirabilis','s__Burkholderia mallei'))

p1= ggplot() +
  geom_jitter( df1,mapping=aes(x=genus, y=hi,fill=type), shape=21,size=0.8 ,stroke=0) +#数据点1
  scale_fill_manual(values=c('#4682B4','#EDEDD9'))+
  geom_segment(df2,mapping=aes(x=genus, xend=genus, y=hic1, yend=hic2), color="black",linewidth=1) +#数据点之间的连线
  geom_point( df2,mapping=aes(x=genus, y=hic1), shape=21,fill='#4682B4',size=4 ,stroke=0.5) +#数据点1
  geom_point( df2,mapping=aes(x=genus, y=hic2),shape=21, fill='#EDEDD9' ,size=4,stroke=0.5) +
  scale_y_continuous(limits=c(10,20),breaks=seq(10,20,2))+
  theme_bw()+theme(axis.text.x =element_text(angle = 45, hjust = 1, vjust = 1),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),legend.position = 'none')+ylab('HI')+xlab(NULL)


p1
p2= ggplot() +
  geom_jitter( df11,mapping=aes(x=species, y=hi,fill=type), shape=21,size=0.8 ,stroke=0) +#数据点1
  scale_fill_manual(values=c('#4682B4','#EDEDD9'))+
  geom_segment(df22,mapping=aes(x=species, xend=species, y=hic1, yend=hic2), color="black",linewidth=1) +#数据点之间的连线
  geom_point( df22,mapping=aes(x=species, y=hic1), shape=21,fill='#4682B4',size=4 ,stroke=0.5) +#数据点1
  geom_point( df22,mapping=aes(x=species, y=hic2),shape=21, fill='#EDEDD9' ,size=4,stroke=0.5) +
  scale_y_continuous(limits=c(10,20),breaks=seq(10,20,2))+
  theme_bw()+theme(axis.text.x =element_text(angle = 45, hjust = 1, vjust = 1),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),legend.position = 'none')+ylab('HI')+xlab(NULL)


p2
p=p2+p1
p
ggsave(file="top10-species-genus-HIc1-HIc2-wide.pdf",device = cairo_pdf, width=6, height=4,p)
