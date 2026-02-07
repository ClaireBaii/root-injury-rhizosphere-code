# Figure 3 分析报告

## 1. 文件清单

| 类型 | 路径 |
|------|------|
| **主脚本** | `figure3_main_alpha_pcoa.R` |
| **补充脚本** | `figure3_distance_to_control.R`（距离到对照），`figure3_rarefaction_robustness.R`（稀释稳健性） |
| **DADA2 脚本** | `run_dada2.R`（16S 测序数据处理） |
| **任务清单** | `ToDO_fig3.md` |
| **输入数据** | `data/完整数据-微生物.csv`（ASV/OTU counts + taxonomy） |
| **输出图文件** | `Fig3_main.pdf / .tiff`，`FigS_distance_to_control.pdf / .tiff`，`FigS_rarefaction_robustness.pdf / .tiff` |
| **统计表格** | `Fig3_alpha_stats.csv`，`Fig3_alpha_dunn_all_pairs.csv`，`Fig3_pcoa_scores.csv`，`Fig3_permanova.csv`，`distance_to_control_stats.csv`，`permanova_robustness.csv` |

---

## 2. Panel 解析

### Panel A — α 多样性（Chao1 / Faith_pd / Shannon / Simpson）

**展示内容**：不同处理组（Control, Treat1–Treat4）的 α 多样性指数比较，包含每个样本点 + 组均值 ± SD，以及与 Control 组的显著性标注。

### Panel B — PCoA（主坐标分析）

**展示内容**：基于 Bray-Curtis 距离的微生物群落 β 多样性分布，显示各处理组在二维空间的聚类情况，并标注 PERMANOVA 统计结果（R²、P 值、置换次数）。

### 补充图 — Distance-to-Control

**展示内容**：量化各处理组到 Control 组的 Bray-Curtis 距离，支持"100% RS 部分回撤"的结论。

### 补充图 — Rarefaction Robustness

**展示内容**：稀释（rarefied）与非稀释（relative abundance）两种处理方式下 PCoA + PERMANOVA 结果的对比，验证结论的稳健性。

---

## 3. 数据来源与样本结构

| 信息项 | 描述 |
|--------|------|
| **数据文件** | `data/完整数据-微生物.csv` |
| **样本命名规则** | `RDYS_{组别编号}-{样本号}`，如 `RDYS_1-1` 为 Control 组第 1 个样本 |
| **组别** | Control（1）, Treat1（2）, Treat2（3）, Treat3（4）, Treat4（5） |
| **样本数** | 每组 n = 3（从样本命名推断；**不确定**：脚本未显式标注样本量） |
| **数据类型** | **样本级别**，每行为一个 ASV/OTU，每列为一个样本的 counts |

---

## 4. 数据预处理流程

### α 多样性

| 步骤 | 方法 |
|------|------|
| 数据读取 | `read.csv(..., fileEncoding = "UTF-8-BOM")` |
| 缺失值处理 | `otu_mat[is.na(otu_mat)] <- 0` |
| 取整 | `round(otu_mat)` |
| 零 ASV 过滤 | 移除所有样本 counts 和为 0 的 ASV |
| Shannon 指数 | 使用 **log2** 底（与坐标轴 ~11-12 一致） |
| Faith's PD | 若无系统发育树，则基于 taxonomy 列构建分类树计算（fallback） |

### β 多样性（PCoA / PERMANOVA）

| 步骤 | 方法 |
|------|------|
| 相对丰度转换 | `sweep(otu_mat, 2, colSums(otu_mat), FUN = "/")` |
| 距离度量 | Bray-Curtis（`vegan::vegdist(..., method = "bray")`） |
| PCoA | `cmdscale(bray, k = 2, eig = TRUE)` |

### Distance-to-Control

| 定义 | 计算方式 |
|------|----------|
| **PairwiseMean** | 样本到所有 Control 样本的 Bray-Curtis 距离均值（Control 自身用 leave-one-out） |
| **Centroid** | 样本到 Control 组相对丰度均值（centroid）的 Bray-Curtis 距离 |

### Rarefaction Robustness

| 处理方式 | 方法 |
|----------|------|
| **Non-rarefied** | 相对丰度 → Bray-Curtis |
| **Rarefied** | `vegan::rrarefy()` 稀释到最小深度 → Bray-Curtis |

---

## 5. 统计/分析方法清单

### 5.1 Kruskal-Wallis 检验

- **目的**：检验各组 α 多样性指数是否存在显著差异
- **函数**：`rstatix::kruskal_test(Value ~ Group)`
- **假设**：非正态分布的连续数据

### 5.2 Dunn 检验（多重比较）

- **目的**：事后两两比较
- **函数**：`rstatix::dunn_test(Value ~ Group, p.adjust.method = "BH")`
- **多重校正**：Benjamini-Hochberg (BH) 法控制 FDR
- **比较方式**：所有组与 Control 组比较

### 5.3 PERMANOVA

- **目的**：检验处理组对微生物群落结构的整体影响
- **函数**：`vegan::adonis2(bray ~ Group, data = meta, permutations = 9999)`
- **距离度量**：Bray-Curtis（相对丰度）
- **置换次数**：9999
- **输出指标**：R²（解释方差比例）、P 值

### 5.4 PCoA（主坐标分析）

- **目的**：降维可视化 β 多样性
- **函数**：`cmdscale(bray, k = 2, eig = TRUE)`
- **解释方差**：PC1/PC2 百分比从特征值计算

---

## 6. 可视化编码要点

### Panel A（α 多样性）

| 元素 | 编码 |
|------|------|
| **散点** | 灰色 jitter 点表示每个样本（`geom_jitter`） |
| **均值** | 彩色圆点表示各组均值（`stat_summary`） |
| **误差线** | 均值 ± 1 SD（`mean_sdl(mult = 1)`） |
| **颜色** | Chao1 红、Faith_pd 绿、Shannon 青、Simpson 紫 |
| **分面** | 按指数分面，Y 轴自由缩放（`facet_wrap(~ Metric, scales = "free_y")`） |
| **显著性** | 仅标注 vs Control 显著的配对（`hide.ns = TRUE`） |
| **标注符号** | `p.adj.signif`（`*`, `**`, `***`, `ns`） |

### Panel B（PCoA）

| 元素 | 编码 |
|------|------|
| **散点** | 各样本位置，颜色按组别映射 |
| **颜色方案** | Control 红、Treat1 绿、Treat2 青、Treat3 蓝、Treat4 紫 |
| **坐标** | 等比例（`coord_equal()`） |
| **副标题** | PERMANOVA 结果（R²、P、置换次数） |

### 补充图（Distance-to-Control / Rarefaction）

| 元素 | 编码 |
|------|------|
| **散点** | 灰色 jitter |
| **均值** | 黑色圆点 |
| **误差线** | 均值 ± 1 SD，黑色 |
| **副标题** | Kruskal-Wallis P 值 |

---

## 7. 关键 R 包/函数

### 统计分析

| 包 | 函数 | 用途 |
|----|------|------|
| **vegan** | `diversity()` | Shannon / Simpson 指数 |
| **vegan** | `vegdist()` | Bray-Curtis 距离矩阵 |
| **vegan** | `adonis2()` | PERMANOVA |
| **vegan** | `rrarefy()` | 稀释采样 |
| **rstatix** | `kruskal_test()` | Kruskal-Wallis 检验 |
| **rstatix** | `dunn_test()` | Dunn 事后检验 |
| **picante** | `pd()` | Faith's PD（需系统发育树） |
| **ape** | `as.phylo()` | 构建系统发育树 |

### 可视化

| 包 | 函数 | 用途 |
|----|------|------|
| **ggplot2** | `geom_jitter()`, `stat_summary()`, `facet_wrap()` | 多面板散点图 |
| **ggpubr** | `stat_pvalue_manual()` | 显著性标注 |
| **patchwork** | `+`, `plot_layout()`, `plot_annotation()` | 图形拼接与标签 |

---

## 8. 潜在风险点

| 风险 | 说明 |
|------|------|
| **样本量不明确** | 脚本从样本命名推断组别，未显式验证每组 n；若实际 n ≠ 预期，可能导致统计功效不足 |
| **Faith's PD 精度** | 若无真实系统发育树，fallback 使用 taxonomy 构建的分类树，精度可能低于 16S 系统发育树 |
| **距离矩阵对称性** | Distance-to-Control 的 leave-one-out 仅影响 Control 组自身，其他组无此处理，可能引入轻微偏差 |
| **稀释深度** | Rarefaction 使用数据集最小深度，若某样本深度极低可能丢失大量信息 |

---

## 9. 输出文件说明

| 文件 | 内容 |
|------|------|
| `Fig3_alpha_stats.csv` | 各指数各组的 n、mean、SD、KW p 值、Dunn q 值 |
| `Fig3_alpha_dunn_all_pairs.csv` | 所有配对的 Dunn 检验结果 |
| `Fig3_pcoa_scores.csv` | 各样本 PC1/PC2 坐标及解释方差 |
| `Fig3_permanova.csv` | PERMANOVA 结果（R²、P、置换次数、距离方法） |
| `distance_to_control_stats.csv` | 两种定义下的均值 ± SD + 检验结果 |
| `permanova_robustness.csv` | 稀释 vs 非稀释 PERMANOVA 对比 |
