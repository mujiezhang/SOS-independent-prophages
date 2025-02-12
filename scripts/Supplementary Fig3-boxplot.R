library(ggplot2)
df<- read.table("Supplementary Fig 3-c-score.txt",header = T)
df$genus<- factor(df$genus, levels = c('g__Dickeya', 'g__Kluyvera', 'g__Yersinia', 'g__Escherichia', 'g__Pectobacterium', 'g__Providencia', 'g__Edwardsiella', 'g__Kosakonia', 'g__Salmonella', 'g__Serratia_J', 'g__Proteus', 'g__Citrobacter', 'g__Erwinia', 'g__Enterobacter', 'g__Morganella', 'g__Serratia', 'g__Cronobacter', 'g__Citrobacter_A', 'g__Leclercia', 'g__Providencia_C', 'g__Klebsiella', 'g__Citrobacter_B', 'g__Pantoea',
                                       'g__Pasteurella','g__Shewanella','g__Vibrio','g__Aeromonas','g__Ralstonia', 'g__Burkholderia', 'g__Cupriavidus','g__Bordetella', 'g__Achromobacter','g__Comamonas',
                                       'g__Pseudomonas', 'g__Pseudomonas_E'))
p=ggplot(df, aes( x = genus,y=score,fill=family))+
  geom_boxplot(outlier.shape = NA)+
 scale_fill_manual(values=c('f__Enterobacteriaceae'='#f5eab0','f__Pasteurellaceae'='#f1c482',
                            'f__Shewanellaceae'='#e8645b','f__Vibrionaceae'='#ed8e67',
                            'f__Aeromonadaceae'='#f6fbcd','f__Burkholderiaceae'='#8472c7',
                            'f__Burkholderiaceae_C'='#7994c0','f__Burkholderiaceae_B'='#a6cbe5'))+
  labs(x = "genus",y = "score") +
  geom_jitter(width=0.1,size=0.5,fill='black')+
  theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())
p

