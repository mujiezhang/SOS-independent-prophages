setwd('E:/课题/4.tem-induce/imeta/simulated_genomes/mock_genome')
#df1=read.table('calculate_rate_distribution.tsv',head=TRUE,sep='\t')
library(ggplot2)
library(dplyr)
library(ggpubr)
library(patchwork)
# 读取数据
data=read.table('calculate_rate_distribution2.tsv',head=TRUE,sep='\t')

hcom_p <- subset(data, type == "h_com" & type2 == "precision")
hcom_f <- subset(data, type == "h_com" & type2 == "success rate")

hcon_p <- subset(data, type == "h_con" & type2 == "precision")
hcon_f <- subset(data, type == "h_con" & type2 == "success rate")

vcom_p <- subset(data, type == "v_com" & type2 == "precision")
vcom_f <- subset(data, type == "v_com" & type2 == "success rate")

vcon_p <- subset(data, type == "v_con" & type2 == "precision")
vcon_f <- subset(data, type == "v_con" & type2 == "success rate")

# 绘图
p1=ggplot(hcom_f, aes(x = percentage, y = ratio)) +
  geom_point(size = 2, color = "steelblue") +  # 散点
  geom_smooth(method = "lm", se = FALSE, color = "firebrick") +  # 线性拟合线
  stat_regline_equation(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~"))) +
  stat_cor(aes(label = paste("P =", after_stat(p))), 
           label.y.npc = 0.1, label.x.npc = 0.1) +
  labs(
    x = "Percentage (%)",
    y = "success rate (%)",
    title = "Host genome completeness"
  ) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 12)
  )

p2=ggplot(hcom_p, aes(x = percentage, y = ratio)) +
  geom_point(size = 2, color = "steelblue") +  # 散点
  geom_smooth(method = "lm", se = FALSE, color = "firebrick") +  # 线性拟合线
  stat_regline_equation(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~"))) +
  stat_cor(aes(label = paste("P =", after_stat(p))), 
           label.y.npc = 0.1, label.x.npc = 0.1) +
  labs(
    x = "Percentage (%)",
    y = "Precision (%)",
    title = "Host genome completeness"
  ) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text( face = "bold", hjust = 0.5))

p3=ggplot(hcon_f, aes(x = percentage, y = ratio)) +
  geom_point(size = 2, color = "steelblue") +  # 散点
  geom_smooth(method = "lm", se = FALSE, color = "firebrick") +  # 线性拟合线
  stat_regline_equation(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~"))) +
  stat_cor(aes(label = paste("P =", after_stat(p))), 
           label.y.npc = 0.1, label.x.npc = 0.1) +
  labs(
    x = "Percentage (%)",
    y = "success rate (%)",
    title = "Host genome contamination"
  ) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text( face = "bold", hjust = 0.5)
  )

p4=ggplot(hcon_p, aes(x = percentage, y = ratio)) +
  geom_point(size = 2, color = "steelblue") +  # 散点
  geom_smooth(method = "lm", se = FALSE, color = "firebrick") +  # 线性拟合线
  stat_regline_equation(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~"))) +
  stat_cor(aes(label = paste("P =", after_stat(p))), 
           label.y.npc = 0.1, label.x.npc = 0.1) +
  labs(
    x = "Percentage (%)",
    y = "Precision (%)",
    title = "Host genome contamination"
  ) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text( face = "bold", hjust = 0.5)
  )

p5=ggplot(vcom_f, aes(x = percentage, y = ratio)) +
  geom_point(size = 2, color = "steelblue") +  # 散点
  geom_smooth(method = "lm", se = FALSE, color = "firebrick") +  # 线性拟合线
  stat_regline_equation(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~"))) +
  stat_cor(aes(label = paste("P =", after_stat(p))), 
           label.y.npc = 0.1, label.x.npc = 0.1) +
  labs(
    x = "Percentage (%)",
    y = "success rate (%)",
    title = "Viral genome completeness"
  ) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5)
  )

p6=ggplot(vcom_p, aes(x = percentage, y = ratio)) +
  geom_point(size = 2, color = "steelblue") +  # 散点
  geom_smooth(method = "lm", se = FALSE, color = "firebrick") +  # 线性拟合线
  stat_regline_equation(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~"))) +
  stat_cor(aes(label = paste("P =", after_stat(p))), 
           label.y.npc = 0.1, label.x.npc = 0.1) +
  labs(
    x = "Percentage (%)",
    y = "Precision (%)",
    title = "Viral genome completeness"
  ) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5)
  )

p7=ggplot(vcon_f, aes(x = percentage, y = ratio)) +
  geom_point(size = 2, color = "steelblue") +  # 散点
  geom_smooth(method = "lm", se = FALSE, color = "firebrick") +  # 线性拟合线
  stat_regline_equation(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~"))) +
  stat_cor(aes(label = paste("P =", after_stat(p))), 
           label.y.npc = 0.1, label.x.npc = 0.1) +
  labs(
    x = "Percentage (%)",
    y = "success rate (%)",
    title = "Viral genome contamination"
  ) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5)
  )

p8=ggplot(vcon_p, aes(x = percentage, y = ratio)) +
  geom_point(size = 2, color = "steelblue") +  # 散点
  geom_smooth(method = "lm", se = FALSE, color = "firebrick") +  # 线性拟合线
  stat_regline_equation(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~"))) +
  stat_cor(aes(label = paste("P =", after_stat(p))), 
           label.y.npc = 0.1, label.x.npc = 0.1) +
  labs(
    x = "Percentage (%)",
    y = "Precision (%)",
    title = "Viral genome contamination"
  ) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5)
  )

layout <- "
AB
CD
EF
GH
"

p=p1 + p2 + p3 +p4 + p5 + p6 +p7 + p8 +  plot_layout(design = layout)
p
ggsave(file="completeness-contamination-with-p-value-20250623.pdf",device = cairo_pdf, width=5, height=10,p)
