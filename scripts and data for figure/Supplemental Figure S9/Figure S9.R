library(ggplot2)

# 读取TSV文件
data0 <- read.table("phage-distri-by-env-percent.tsv",head=TRUE,sep='\t')

# 确保环境因子顺序（按原始顺序）
data0$env <- factor(data0$env, levels = c('human','animal','plant','soil','rhizosphere','sediment','freshwater','wastewater','seawater','other'))
data0$type <- factor(data0$type, levels = c('SiPs','SuPs','SdPs'))

# 创建自定义颜色方案（可选）
type_colors <- c(SiPs = "#E3B680", SuPs = "#DBDCDC", SdPs = "#B2CCE5")

# 绘制百分比堆积柱状图
p=ggplot(data, aes(x = env, y = num, fill = type)) +
  geom_col(position = "stack", width = 0.8) +
  scale_fill_manual(values = type_colors) +
  labs(x = NULL,
       y = "Percentage (%)") +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "top"
  ) +
  geom_text(aes(label = sprintf("%.1f", num)), 
            position = position_stack(vjust = 0.5),
            size = 3, color = "black")
p

ggsave(file="env-distribution.pdf",device = cairo_pdf, width=4, height=5,p)

data <- read.table("human-distri-by-env-percent.tsv",head=TRUE,sep='\t')
data2 <- read.table("animal-distri-by-env-percent.tsv",head=TRUE,sep='\t')
# 确保环境因子顺序（按原始顺序）
data$env <- factor(data$env, levels = c('stool','blood','urine','sputum','gut','wound','other'))
data$type <- factor(data$type, levels = c('SiPs','SuPs','SdPs'))
data2$env <- factor(data2$env, levels = c('bovine','poultry','pig','fish','dog','bird','murine','other'))
data2$type <- factor(data2$type, levels = c('SiPs','SuPs','SdPs'))
# 创建自定义颜色方案（可选）
type_colors <- c(SiPs = "#E3B680", SuPs = "#DBDCDC", SdPs = "#B2CCE5")

# 绘制百分比堆积柱状图
p1=ggplot(data, aes(x = env, y = num, fill = type)) +
  geom_col(position = "stack", width = 0.8) +
  scale_fill_manual(values = type_colors) +
  labs(x = NULL,
       y = "Percentage (%)") +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "top"
  ) +
  geom_text(aes(label = sprintf("%.1f", num)), 
            position = position_stack(vjust = 0.5),
            size = 3, color = "black")
p1
p2=ggplot(data2, aes(x = env, y = num, fill = type)) +
  geom_col(position = "stack", width = 0.8) +
  scale_fill_manual(values = type_colors) +
  labs(x = NULL,
       y = "Percentage (%)") +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "top"
  ) +
  geom_text(aes(label = sprintf("%.1f", num)), 
            position = position_stack(vjust = 0.5),
            size = 3, color = "black")
p2
ggsave(file="human-distribution.pdf",device = cairo_pdf, width=4, height=5,p1)
ggsave(file="animal-distribution.pdf",device = cairo_pdf, width=4, height=5,p2)
