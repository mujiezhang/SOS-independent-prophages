# 加载所需的库
library(ggplot2)
library(readr)
library(ggsignif)
library(gridExtra)  # 用于图形排布
library(grid)       # 用于图形处理
library(gtable)     # 用于图形表格操作

# 辅助函数：提取 ggplot 对象的图例
get_legend <- function(myplot) {
  tmp <- ggplotGrob(myplot)
  leg <- gtable::gtable_filter(tmp, "guide-box")
  return(leg)
}

# 定义一个函数来生成单个图形，并保留所有注释
create_violin_plot <- function(file, title = "Coding Density", y_label = "Coding density (kb⁻¹)", show_legend = FALSE) {
  
  # 读取合并后的文件
  combined_data <- read_delim(file, delim = "\t")
  
  # 检查读取的数据
  print("数据预览:")
  print(head(combined_data))

  names(combined_data) <- c("name", "value", "Source")
  
  # 将 Source 转换为因子，并设置顺序（仅包含两个类别）
  combined_data <- combined_data[combined_data$Source %in% c("SiPs", "SdPs"), ]
  combined_data$Source <- factor(combined_data$Source, levels = c("SiPs", "SdPs"))
  
  # 检查合并后的数据框
  print("合并后的数据预览:")
  print(head(combined_data))
  
  # 计算两两之间的 p 值（双侧ttest检验）
  p_values <- list()
  # 仅计算 "SdPs" 和 "SiPs" 之间的 p 值
  p_values[["SdPs_vs_SiPs"]] <- t.test(value ~ Source, data = subset(combined_data, Source %in% c("SdPs", "SiPs")), alternative = "two.sided")$p.value
  print("原始 p 值:")
  print(p_values)
  
  
  # 创建一个函数来将 p 值转换为显著性符号
  get_significance <- function(p) {
    if (is.na(p)) {
      return("")
    }
    
    if (p < 0.001) {
      return("***")
    } else if (p < 0.01) {
      return("**")
    } else if (p < 0.05) {
      return("*")
    } else {
      return("")
    }
  }
  
  # 将调整后的 p 值转换为显著性符号
  significance <- sapply(p_values, get_significance)
  
  # 自定义颜色（仅包含两个类别）
  custom_colors <- c(
    "SdPs" = "#accde8", 
    "SiPs" = "#ECB477"
  )
  
  # 动态计算 y_position 以适应数据范围（仅一个比较）
  y_max <- max(combined_data$value, na.rm = TRUE)
  y_positions <- y_max * 1.05 
  
  # 绘制小提琴图并显示散点及显著性标记
  p <- ggplot(combined_data, aes(x = Source, y = value, fill = Source)) +  
    geom_violin(alpha = 1, width = 0.6, trim = TRUE, draw_quantiles = c(0.25, 0.5, 0.75)) +  
    scale_fill_manual(values = custom_colors) +  
    labs(title = title, x = "Source", y = y_label) +  
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),  
      panel.border = element_rect(color = "black", fill = NA, size = 0.5), 
      axis.text.x = element_text(angle = 0, hjust = 1), 
      axis.ticks.length.y = unit(0.05, "cm"), 
      axis.ticks.y = element_line(color = "black"),
      legend.position = ifelse(show_legend, "right", "none") 
    ) +
    
    # 使用 ggsignif 包添加显著性标记
    geom_signif(
      comparisons = list(
        c("SdPs", "SiPs") 
      ),
      map_signif_level = TRUE,
      textsize = 3,
      color = "black",
      y_position = y_positions, 
      tip_length = 0, 
      test = "t.test",  
      annotations = significance,  
      vjust = 0.6
    ) +
    geom_text(
      aes(
        x = 1.5, 
        y = y_positions * 1.08, 
        label = paste0("p=", signif(p_values[["SdPs_vs_SiPs"]], 3))  
      ),
      color = "red", size = 2
    ) +
    # 添加平均值标记（三角形）和平均值标签
    stat_summary(fun = mean, geom = "point", shape = 17, size = 1, color = "#e05c00", fill = "#e05c00") + 
    stat_summary(fun = mean, geom = "text", aes(label = round(after_stat(y), 3)), vjust = -0.5, color = "#e05c00") +  
    
    # 使用 geom_line() 添加灰色阴影线
    stat_summary(
      fun = mean, 
      geom = "line", 
      aes(group = 1),  
      color = "#9cc3d4", 
      linetype = "solid",
      size = 3,  
      alpha = 0.4  
    ) +
    
    # 使用 geom_line() 添加虚线
    stat_summary(
      fun = mean, 
      geom = "line", 
      aes(group = 1),  
      color = "#c00000",  
      linetype = "dashed", 
      size = 0.5 
    )

  return(p)
}

# 定义多个数据集的文件路径以及对应的标题和纵轴标签
data_sets <- list(
  list(
    file = "3_all_d2_k6.txt", 
    title = "Dataset 1: d2",
    y_label = "d2"
  ),
  list(
    file = "4_all_d2shepp_k6.txt",
    title = "Dataset 2: d2shepp",
    y_label = "d2shepp"
  ),
  list(
    file = "5_all_Eu_k6.txt", 
    title = "Dataset 3: Eu",
    y_label = "Eu"
  ),
  list(
    file = "6_all_EuF_k6.txt",  
    title = "Dataset 4: EuF",
    y_label = "EuF"
  ),
  list(
    file = "7_all_Hao_k6.txt", 
    title = "Dataset 5: Hao",
    y_label = "Hao"
  ),
  list(
    file = "8_all_JS_k6.txt", 
    title = "Dataset 6: JS",
    y_label = "JS"
  ),
  list(
    file = "9_all_Ma_k6.txt", 
    title = "Dataset 7: Ma",
    y_label = "Ma"
  ),
  list(
    file = "10_all_Ch_k6.txt",  
    title = "Dataset 8: Ch",
    y_label = "Ch"
  )
)

# 检查 data_sets 列表是否正确
print("定义的数据集列表:")
print(data_sets)

# 生成所有图形并存储在列表中
plot_list <- lapply(seq_along(data_sets), function(i) {
  ds <- data_sets[[i]]
  # 仅对第一个图形显示图例，其余图形隐藏图例
  create_violin_plot(
    file = ds$file, 
    title = ds$title, 
    y_label = ds$y_label,
    show_legend = ifelse(i == 1, TRUE, FALSE)
  )
})

# 提取第一个图形的图例
legend <- get_legend(plot_list[[1]])

# 移除所有图形的图例（已经在生成图形时通过参数控制）
plot_list_no_legend <- lapply(plot_list, function(p) p + theme(legend.position = "none"))

# 使用 grid.arrange 排布图形和图例
grid.arrange(
  arrangeGrob(grobs = plot_list_no_legend, ncol = 4, nrows = 3), 
  legend, 
  #ncol = 5, 
  widths = c(20, 1) ,
  heights = c(4, 1) 
)

