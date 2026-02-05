### 📉 图4：关键分泌物 log2FC 趋势图（Top50）+ 回归线 + 显著性标注（线性趋势）
# 文件：top50_log2fc_trend_data.csv
# 包依赖：ggplot2, dplyr, broom

library(ggplot2)
library(dplyr)
library(broom)
setwd("E:/我的/蔡金秀论文/榆树根系/断根组数据")
trend_df <- read.csv("top50_log2fc_trend_data.csv")

# 计算每个分泌物的线性拟合显著性（p值）
trend_df$Treatment <- factor(trend_df$Treatment,
                             levels = c("Control", "Processing_1", "Processing_2", "Processing_3", "Processing_4"))
trend_df$Treatment_numeric <- as.numeric(trend_df$Treatment)

p_table <- trend_df %>%
  group_by(Compound) %>%
  filter(!is.na(log2FC)) %>%
  filter(n() >= 3) %>%
  do({
    fit <- tryCatch(lm(log2FC ~ Treatment_numeric, data = .), error = function(e) NULL)
    if (!is.null(fit)) tidy(fit) else tibble(term = NA, p.value = NA)
  }) %>%
  filter(term == "Treatment_numeric") %>%
  select(Compound, p.value)

# 标记显著性（p < 0.05）

# 导出按 p 值排序的显著性表
p_table_sorted <- p_table %>% arrange(p.value)
write.csv(p_table_sorted, "trend_linear_significance_table.csv", row.names = FALSE)
trend_df <- left_join(trend_df, p_table, by = "Compound")
trend_df$p_label <- ifelse(trend_df$p.value < 0.05, "*", "")

trend_df$Treatment <- factor(trend_df$Treatment,
                             levels = c("Control", "Processing_1", "Processing_2", "Processing_3", "Processing_4"))

ggplot(trend_df, aes(x = Treatment, y = log2FC, group = Compound, color = Compound)) +
  geom_point(size = 2) +
  geom_line(size = 1) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", size = 0.8) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none") +
  labs(title = "Log2 Fold Change Trend of Top 50 Root Exudates",
       y = "log2FC vs Control",
       x = "Treatment Group") +
  facet_wrap(~Compound, scales = "free_y", ncol = 5) +
  geom_text(data = subset(trend_df, Treatment == "Processing_4"),
            aes(x = Treatment, y = max(log2FC, na.rm = TRUE), label = p_label),
            inherit.aes = FALSE, size = 3, vjust = -0.5)