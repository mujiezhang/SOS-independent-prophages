library(ggplot2)
library(patchwork)

# 读取TSV文件
data <- read.table("retained-top-20-for-r.tsv",head=TRUE,sep='\t')
data2 <- read.table("filtered-out-top-20-for-r.tsv",head=TRUE,sep='\t')
# 按原始顺序
data$genus <- factor(data$genus, levels = unique(data$genus))
data$reason <- factor(data$reason, levels = c('no prophage','unqualified prophage','no lexa','no box','fail','keep'))
data2$genus <- factor(data2$genus, levels = unique(data2$genus))
data2$reason <- factor(data2$reason, levels = c('no prophage','unqualified prophage','no lexa','no box','fail'))

# 创建自定义颜色方案（可选）
type_colors <- c('no prophage' = "#E3B680", 'unqualified prophage' = "#DBDCDC", 'no lexa' = "#CBE3A6",
                 'no box' = "#949B93", 'fail' = "#93B9BE", 'keep' = "#B2CCE5")

# 绘制百分比堆积柱状图
p1=ggplot(data, aes(x = genus, y = num, fill = reason)) +
  geom_col(position = "fill", width = 0.8) +
  scale_fill_manual(values = type_colors) +
  labs(x = NULL,
       y = "Percentage (%)") +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "top"
  ) 

p1
p2=ggplot(data2, aes(x = genus, y = num, fill = reason)) +
  geom_col(position = "fill", width = 0.8) +
  scale_fill_manual(values = type_colors) +
  labs(x = NULL,
       y = "Percentage (%)") +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "top"
  ) 
p2
p=p1/p2
p
ggsave(file="retained-and-filtered-out.pdf",device = cairo_pdf, width=5, height=6,p)
