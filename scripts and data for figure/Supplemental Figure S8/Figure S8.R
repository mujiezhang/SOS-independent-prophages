library(ggplot2)
library(ggsignif)
df<- read.table("sdp_sip.tsv",header = T)
df$type <- factor(df$type, levels = c('sip','sdpvssip',"sdp"))
mycomparisons1 <- list(c("sdp",'sip'),c("sdp", "sdpvssip"),c("sdpvssip",'sip'))
p=ggplot(df, aes( x = type,y=wgrr,fill = type))+
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
  geom_signif(
    comparisons = mycomparisons1,
    test =t.test,
    y_position = c(105,110,115),
    tip_length = 0.01,
  )+
  geom_point(data = data.frame(type = c('sdp', 'sdpvssip', 'sip'), wgrr = c(8.65, 7.79, 12.86)),
             aes(x = type, y = wgrr), color = "red", size = 3)+
  scale_fill_manual(values=c('sdp'='#ACCDE8','sip'='#ECB477','sdpvssip'='#ECF4DD'))+
  labs(x = NULL,
       y = "Percentage") +
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),legend.position ='none')
p
ggsave(file="violin-wgrr.pdf",device = cairo_pdf, width=3, height=5,p)
