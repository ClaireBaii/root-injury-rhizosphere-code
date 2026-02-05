# 重新读取数据（确保清洁）
setwd("E:/我的/蔡金秀论文/榆树根系/断根组数据")
phylum_df <- read.csv("断根 phylum(1).csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)

# 必须预先执行过：phylum_df 已加载并包含 Ratio 列
# 若未加载，请添加：
# phylum_df <- read.csv("断根 phylum(1).csv", header = TRUE, row.names = 1)

library(dplyr)
library(tidyr)
library(ggplot2)

# 添加截断比变量（如尚未加）
phylum_df$Ratio <- c(0, 25, 50, 75, 100)

# 1️⃣ Spearman ρ 分析
rho_table <- sapply(phylum_df[, !colnames(phylum_df) %in% c("Group", "Ratio")], function(x) {
  cor.test(x, phylum_df$Ratio, method = "spearman")$estimate
})
pval_table <- sapply(phylum_df[, !colnames(phylum_df) %in% c("Group", "Ratio")], function(x) {
  cor.test(x, phylum_df$Ratio, method = "spearman")$p.value
})

spearman_result <- data.frame(Phylum = names(rho_table),
                              Rho = as.numeric(rho_table),
                              p_value = as.numeric(pval_table))

# 2️⃣ 筛选趋势型响应菌门（|ρ| > 0.7）
trend_phyla <- spearman_result %>% filter(abs(Rho) > 0.7)

# 3️⃣ 清洗函数（统一用于列名和菌门名）
clean_colnames <- function(x) {
  x <- trimws(x)
  x <- gsub("[\r\n]", "", x)
  x <- gsub("　", "", x)
  x <- gsub("\\.rho$", "", x)  # 处理 .rho 后缀
  x <- gsub("\\.+$", "", x)    # 去除末尾多余句点
  return(x)
}

# 4️⃣ 应用清洗：列名、菌门名统一处理
colnames(phylum_df) <- clean_colnames(colnames(phylum_df))
trend_phyla$Phylum <- clean_colnames(trend_phyla$Phylum)

# 5️⃣ 匹配趋势菌门列名并绘图
if (nrow(trend_phyla) == 0) {
  cat("🚫 无趋势性响应菌门（|ρ| > 0.7），跳过绘图。\n")
} else {
  real_phylum_cols <- setdiff(colnames(phylum_df), c("Ratio", "Group"))
  trend_cols <- intersect(trend_phyla$Phylum, real_phylum_cols)
  
  missing_cols <- setdiff(trend_phyla$Phylum, colnames(phylum_df))
  if (length(trend_cols) == 0) {
    cat("⚠️ 无法匹配任何趋势菌门名称，请检查列名是否异常。\n")
    print("未匹配成功的菌门名如下：")
    print(missing_cols)
  } else {
    # 构建绘图数据
    plot_df <- phylum_df[, trend_cols, drop = FALSE]
    plot_df$Sample <- rownames(phylum_df)
    plot_df$Group <- c("Control", "Treat25", "Treat50", "Treat75", "Treat100")
    
    # 转为长表
    plot_df <- plot_df %>%
      pivot_longer(-c(Sample, Group), names_to = "Phylum", values_to = "Abundance")
    
    # 绘图
    p <- ggplot(plot_df, aes(x = Group, y = Abundance, fill = Group)) +
      geom_bar(stat = "identity", position = "dodge") +
      facet_wrap(~ Phylum, scales = "free_y") +
      scale_fill_brewer(palette = "Set2") +
      labs(x = "主根截断组", y = "相对丰度", title = "趋势响应型菌门（|ρ| > 0.7）在各组丰度分布") +
      theme_minimal(base_size = 13) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            strip.text = element_text(face = "bold"))
    
    print(p)
  }
}
# 假设分泌物表已加载为 exudate_df，行为样本，列为代谢物
exudate_df$Ratio <- c(0, 25, 50, 75, 100)
colnames(exudate_df)[1:10]
# 以某一类分泌物（如 Fatty Acyls）为例绘图
library(tidyr)
library(dplyr)
ggplot(exudate_df, aes(x = Ratio, y = `Fatty.Acyls`)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 1) +
  theme_minimal(base_size = 13) +
  labs(x = "主根截断比例 (%)", y = "Fatty Acyls 丰度",
       title = "Fatty Acyls 随主根截断比例的变化趋势")


# 构建长表格（提取你感兴趣的分泌物列）
df_plot <- exudate_df %>%
  select(Ratio, Fatty.Acyls, Phenols, Tannins = Diazines, Alkaloids = Heteroaromatic.compounds) %>%
  pivot_longer(-Ratio, names_to = "Metabolite", values_to = "Abundance")

# 绘图
ggplot(df_plot, aes(x = Ratio, y = Abundance, color = Metabolite)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal(base_size = 13) +
  labs(x = "主根截断比 (%)", y = "丰度", title = "主要分泌物随截断比的变化趋势")
