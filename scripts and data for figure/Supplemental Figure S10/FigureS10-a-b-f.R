library(ggplot2)
library(patchwork)
df1=read.table('sdp-sip.txt',head=TRUE,sep='\t')
df2=read.table('key-gene-data.tsv',head=TRUE,sep='\t')
df2$type <- factor(df2$type,levels =  c('SiPs-L','SdPs-L','SiPs-S','SdPs-S'))
df2$functio <- factor(df2$functio,levels = c('integrase','terminase','portal-associated','lysis-associated','tail- and baseplate-associated'))
df1$type2 <- factor(df1$type2,levels = c('SiPs','SdPs'))
df1$type <- factor(df1$type,levels =  c('SiPs-L','SdPs-L','SiPs-S','SdPs-S'))
p<-ggplot(df1, aes(x = contig_length,fill=type2)) +
  scale_x_continuous(limits=c( 0 , 100000 ),breaks = seq(0, 100000, by = 15000))+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),legend.position= 'none')+
  geom_density(alpha=0.5)+ xlab('Virus length (bp)')+ylab('Density')+
  scale_fill_manual(values=c("#ECB477","#accde8"))+
  facet_grid(type2 ~ .) 
p
p2=ggplot(df1, aes( x = type,y=completeness,fill=type))+
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values=c("#ECB477","#accde8","#ead1b9","#dfe8ef"))+
  labs(x = NULL,
       y = "Completeness (%)") +
  geom_jitter(width=0.1,size=0.2,fill='black',alpha=0.2,stroke=0)+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),legend.position = 'none')
p2
p3=ggplot(df2,aes(x=functio,y=percent,fill=type))+
  geom_bar(stat = 'identity', 
           position = 'dodge', 
           width = 0.8, 
           color='black')+xlab('Protein function')+ylab('Percentage (%)')+ 
  scale_fill_manual(values=c("#ECB477","#accde8","#ead1b9","#dfe8ef"))+
  geom_text(aes(label=percent),size=3,position = position_dodge(width = 0.8),vjust=-0.3)+   
  theme_bw()+theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()) 
p3
p1=p/p2/p3
p1
ggsave(file="length+completeness+protein-function.pdf",device = cairo_pdf, width=4, height=5,p1)

