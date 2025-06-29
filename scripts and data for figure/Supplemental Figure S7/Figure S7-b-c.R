library(ggplot2)
library(ggsignif)

df=read.table('multi-phage-for-r.tsv',head=TRUE,sep='\t')

df$type<- factor(df$type, levels = c('2','3','4','5',
                                     '6','7','8','9',
                                     '10','11','12','13',
                                     '14','15','16','17','18','19'))

p=ggplot(df, aes( x = type,y=num))+
  geom_col(width = 0.7,colour='black',linewidth=0.2,fill='#adcde8')+
  ylim(0,3500)+ylab('Genomes')+xlab('Phages')+
  geom_text(aes(label = num), vjust = -1)+
  theme_bw()+theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = 'none')
p

ggsave(file="multi-phage.pdf",device = cairo_pdf, width=6, height=3,p)

df1<- read.table("regulation-consistence-for-r.tsv",header = T)

df1$genus<- factor(df1$genus, levels = c('g__Dickeya', 'g__Kluyvera', 'g__Yersinia', 'g__Escherichia', 
                                         'g__Pectobacterium', 'g__Providencia', 'g__Edwardsiella', 
                                         'g__Kosakonia', 'g__Salmonella', 'g__Serratia_J', 'g__Proteus',
                                         'g__Citrobacter', 'g__Erwinia', 'g__Enterobacter', 'g__Morganella', 
                                         'g__Serratia', 'g__Cronobacter', 'g__Citrobacter_A', 'g__Leclercia',
                                         'g__Providencia_C', 'g__Klebsiella', 'g__Citrobacter_B', 'g__Pantoea',
                                       'g__Pasteurella','g__Shewanella','g__Vibrio','g__Aeromonas',
                                       'g__Ralstonia', 'g__Burkholderia', 'g__Cupriavidus','g__Bordetella', 'g__Achromobacter','g__Comamonas',
                                       'g__Pseudomonas', 'g__Pseudomonas_E'))

p1=ggplot(df1, aes( x = genus,y=score,fill=family))+
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values=c('f__Enterobacteriaceae'='#f5eab0','f__Pasteurellaceae'='#f1c482',
                             'f__Shewanellaceae'='#e8645b','f__Vibrionaceae'='#ed8e67',
                             'f__Aeromonadaceae'='#f6fbcd','f__Burkholderiaceae'='#8472c7',
                             'f__Burkholderiaceae_C'='#7994c0','f__Burkholderiaceae_B'='#a6cbe5'))+
  labs(x = "genus",
       y = "score") +
  geom_jitter(width=0.1,size=0.5,fill='black')+
  theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())
p1
ggsave(file="boxplot.pdf",device = cairo_pdf, width=15, height=5,p1)
