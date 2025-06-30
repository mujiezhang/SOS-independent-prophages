# 加载所需的库
library(ggplot2)
library(dplyr)
library(gridExtra)  # 确保已安装，用于 grid.arrange

# 读取数据
data <- read.delim("6_growth.tsv", sep = "\t")

# 计算百分比
data <- data %>%
  group_by(Category) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  ungroup()  # 取消分组，避免后续操作的潜在问题

# 设置 Category 的顺序
data$Category <- factor(data$Category, levels = c("SiTVs", "SuTVs", "SdTVs"))

# 设置 Group 的顺序，确保 "slow" 在 "fast" 之后
data$Group <- factor(data$Group, levels = c("slow", "fast"))

# 为 fast 组和 slow 组分别设置颜色
data <- data %>%
  mutate(Group_Color = ifelse(Group == "fast", as.character(Category), "slow"))

# 将 Group_Color 转换为因子，并确保 'slow' 是最后一个水平
data$Group_Color <- factor(data$Group_Color, levels = c("slow", "SiTVs", "SuTVs", "SdTVs"))

# 定义颜色映射
color_mapping <- c(
  "SiTVs" = "#ECB477",
  "SuTVs" = "#ECF4DD", 
  "SdTVs" = "#accde8", 
  "slow" = "gray"      
)

# 绘制百分比堆积条形图
p <- ggplot(data, aes(x = Category, y = Percentage, fill = Group_Color)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Growth Rate Distribution by Category",
    x = "Category",
    y = "Percentage (%)",
    fill = "Growth Rate"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.ticks.y = element_line(color = "black"),  
    axis.line.y = element_line(color = "black")  
  ) +
  scale_fill_manual(values = color_mapping) 

print(p)

