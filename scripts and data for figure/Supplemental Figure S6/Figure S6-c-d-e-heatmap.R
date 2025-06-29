df1 <- read.table('host-HI-infor-top10genus-with-family.tsv', header = TRUE, sep = '\t')
df2<- read.table('host-HI-infor-top10species-with-genus.tsv', header = TRUE, sep = '\t')
library(ggplot2)
library(reshape2)
library(patchwork)
# 设置genus因子顺序
species_order <- c('s__Escherichia coli', 's__Klebsiella pneumoniae',
                   's__Klebsiella variicola','s__Klebsiella quasipneumoniae','s__Salmonella enterica',
                   's__Pseudomonas aeruginosa','s__Enterobacter hormaechei_A',
                   's__Citrobacter freundii','s__Proteus mirabilis','s__Burkholderia mallei')
genus_order <- c('g__Escherichia','g__Klebsiella','g__Salmonella','g__Enterobacter','g__Citrobacter','g__Pseudomonas',
                 'g__Pseudomonas_E','g__Vibrio','g__Burkholderia','g__Aeromonas')

# 设置family因子顺序
family_order <- c('f__Enterobacteriaceae', 'f__Pseudomonadaceae', 
                  'f__Vibrionaceae', 'f__Burkholderiaceae', 'f__Aeromonadaceae')

# 应用排序
df1$genus <- factor(df1$genus, levels = genus_order)
df1$family <- factor(df1$family, levels = family_order)
df2$species <- factor(df2$species, levels = species_order)
# 函数：创建简洁的p值热图

create_clean_pheatmap <- function(data, group_var) {
  # 分别提取hic1和hic2数据
  hic1_data <- subset(data, type == "hic1")
  hic2_data <- subset(data, type == "hic2")
  
  # 获取分组顺序
  group_levels <- levels(data[[group_var]])
  n <- length(group_levels)
  
  # 创建空矩阵存储p值
  hic1_pvalues <- matrix(NA, nrow = n, ncol = n)
  hic2_pvalues <- matrix(NA, nrow = n, ncol = n)
  rownames(hic1_pvalues) <- colnames(hic1_pvalues) <- group_levels
  rownames(hic2_pvalues) <- colnames(hic2_pvalues) <- group_levels
  
  # 使用Welch t检验计算所有两两比较的p值
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      group1 <- group_levels[i]
      group2 <- group_levels[j]
      
      # 对hic1数据执行Welch t检验
      data1_i <- hic1_data$hi[hic1_data[[group_var]] == group1]
      data1_j <- hic1_data$hi[hic1_data[[group_var]] == group2]
      if (length(data1_i) > 1 && length(data1_j) > 1) {
        test1 <- t.test(data1_i, data1_j, var.equal = FALSE)  # Welch t检验
        hic1_pvalues[j, i] <- test1$p.value
      }
      
      # 对hic2数据执行Welch t检验
      data2_i <- hic2_data$hi[hic2_data[[group_var]] == group1]
      data2_j <- hic2_data$hi[hic2_data[[group_var]] == group2]
      if (length(data2_i) > 1 && length(data2_j) > 1) {
        test2 <- t.test(data2_i, data2_j, var.equal = FALSE)  # Welch t检验
        hic2_pvalues[i, j] <- test2$p.value
      }
    }
  }
  
  # 创建组合矩阵（下三角=hic1，上三角=hic2）
  combined_matrix <- matrix(NA, nrow = n, ncol = n)
  rownames(combined_matrix) <- colnames(combined_matrix) <- group_levels
  
  # 填充下三角（hic1数据）
  combined_matrix[lower.tri(combined_matrix, diag = FALSE)] <- 
    hic1_pvalues[lower.tri(hic1_pvalues, diag = FALSE)]
  
  # 填充上三角（hic2数据）
  combined_matrix[upper.tri(combined_matrix, diag = FALSE)] <- 
    hic2_pvalues[upper.tri(hic2_pvalues, diag = FALSE)]
  
  # 转换为长格式数据框
  p_df <- melt(combined_matrix, na.rm = TRUE)
  names(p_df) <- c("Group1", "Group2", "P_Value")
  
  # 创建显著性标签
  p_df$Significance <- ifelse(p_df$P_Value < 0.001, "***",
                              ifelse(p_df$P_Value < 0.01, "**",
                                     ifelse(p_df$P_Value < 0.05, "*", "")))
  
  # 确定颜色标度的最大值
  max_p <- max(p_df$P_Value, na.rm = TRUE)
  if(max_p >= 0.1) {
    max_p <- ceiling(max_p * 10) / 10
  } else if (max_p < 0.1) {
    max_p <- 0.1
  }
  
  # 创建热图
  p_heatmap <- ggplot(p_df, aes(x = Group2, y = Group1, fill = P_Value)) +
    geom_tile(color = "white", size = 0.8) +
    geom_text(aes(label = Significance), size = 6, color = "black", fontface = "bold") +
    scale_fill_gradientn(
      colors = c("#F46D43", "white", "#00B7FF"),
      values = scales::rescale(c(0, 0.05, max_p)),
      limits = c(0, max_p),
      na.value = "grey90"
    ) +
    scale_x_discrete(position = "top", limits = group_levels) +
    scale_y_discrete(limits = rev(group_levels)) +
    theme_minimal() +
    theme(
      axis.text.x.top = element_text(angle = 45, hjust = 0, size = 10, color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      panel.border = element_rect(color = "gray", fill = NA, size = 0.5),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    ) +
    coord_fixed() +
    labs(fill = "P-Value") +
    guides(fill = guide_colorbar(barheight = unit(2.5, "cm")))
  
  return(p_heatmap)
}
# 创建genus水平热图
s <- create_clean_pheatmap(df2, "species") +
  ggtitle("Species Level Comparisons")

g <- create_clean_pheatmap(df1, "genus") +
  ggtitle("Genus Level Comparisons")

# 创建family水平热图
f <- create_clean_pheatmap(df1, "family") +
  ggtitle("Family Level Comparisons")

# 显示并保存图形
print(s)
print(g)
print(f)
p=s/g/f
print(p)
ggsave("species_genue-family-p-value-welch.pdf", plot = p, 
       device = cairo_pdf, width = 6, height = 15)
