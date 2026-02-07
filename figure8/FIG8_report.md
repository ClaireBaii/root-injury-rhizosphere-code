# Figure 8 报告：Procrustes 分析（代谢组 vs 微生物组）

---

## 1. 文件清单

| 类型 | 文件名 | 路径 |
|------|--------|------|
| **主脚本** | `figure8 Procrustes 分析与 PROTEST 检验（美化）.R` | `./figure8/` |
| **ToDO 文档** | `ToDO_fig8.md` | `./figure8/` |
| **输入数据 1** | `完整数据-分泌物.csv` | `./data/` |
| **输入数据 2** | `完整数据-微生物.csv` | `./data/` |
| **输出图 (主图)** | `Fig8_combined.pdf`, `Fig8_combined.tiff` | `./figure8/` |
| **输出图 (备用)** | `Fig8_main.pdf`, `FigS8_procrustes_residuals.pdf` | `./figure8/` |

---

## 2. 图的主旨与 Panel 解释

### Panel A：Procrustes Alignment

**展示内容**：代谢组 (Exudates) 与微生物组 (Microbiota) 两组数据在降维空间中的 Procrustes 对齐结果。

**回答的问题**：根系分泌物代谢谱与根际微生物群落组成之间是否存在显著的协同变化模式（即两者的群落结构是否匹配）？

**关键元素**：
- 圆形点 (●)：代谢组 PCA 坐标（目标矩阵）
- 三角形点 (▲)：微生物组 PCoA 坐标（旋转后拟合）
- 箭头：连接同一样本的两种数据类型，箭头长度 = Procrustes 残差
- 颜色编码：不同处理比例（0%, 25%, 50%, 75%, 100%）

### Panel B：Procrustes Residuals

**展示内容**：各处理组的 Procrustes 残差分布（箭头长度）。

**回答的问题**：不同处理组中，代谢组与微生物组的匹配程度是否存在差异？残差越小表示两组数据协同性越好。

---

## 3. 数据来源与样本结构

### 3.1 代谢组数据（`完整数据-分泌物.csv`）

| 属性 | 值 |
|------|-----|
| 样本命名 | `YS[比例]_F-[重复]`，如 `YS0_F-1`, `YS100_F-3` |
| 处理组 | 0%, 25%, 50%, 75%, 100%（共 5 组） |
| 每组重复 | 3 个（基于命名推断） |
| 总样本数 | 15（5 组 × 3 重复） |
| 数据类型 | 代谢物丰度矩阵（行 = 化合物，列 = 样本） |

### 3.2 微生物组数据（`完整数据-微生物.csv`）

| 属性 | 值 |
|------|-----|
| 样本命名 | `RDYS_[组号]-[重复]`，如 `RDYS_1-1`, `RDYS_5-3` |
| 命名映射 | RDYS_1 → YS0%, RDYS_2 → YS25%, ..., RDYS_5 → YS100% |
| 总样本数 | 15（与代谢组配对） |
| 数据类型 | OTU/ASV counts + taxonomy |

### 3.3 样本配对

脚本通过 `map_rdys_to_ys()` 函数将微生物组样本名映射到代谢组样本名，确保两组数据的样本一一对应。最终取交集样本进行 Procrustes 分析。

---

## 4. 数据预处理流程

### 4.1 代谢组预处理

```r
# 1. 读取数据并提取样本列（正则匹配 YS\d+_F-\d+）
# 2. 转置矩阵（样本 × 化合物）
# 3. log2(x + 1) 转换
mat <- log2(mat + 1)

# 4. 过滤零方差特征（sd = 0 或 NA）
sds <- apply(mat, 2, sd)
keep <- is.finite(sds) & (sds > 0)
mat <- mat[, keep]

# 5. PCA（中心化 + 标准化）
pca <- prcomp(mat, center = TRUE, scale. = TRUE)
scores <- pca$x[, 1:2]  # 提取前2个主成分
```

### 4.2 微生物组预处理

```r
# 1. 读取数据并提取样本列（正则匹配 RDYS_\d+-\d+）
# 2. 过滤全零行（无检出的 OTU）
counts <- counts[rowSums(counts) > 0, ]

# 3. 转置矩阵（样本 × OTU）
# 4. 相对丰度标准化（TSS，按样本总和归一化）
rel <- sweep(counts_t, 1, rowSums(counts_t), "/")

# 5. Bray-Curtis 距离矩阵
bray <- vegdist(rel, method = "bray")

# 6. PCoA（经典多维标度法）
pcoa <- cmdscale(bray, k = 2, eig = TRUE)
# 如果有负特征值问题，添加校正：
pcoa <- cmdscale(bray, k = 2, eig = TRUE, add = TRUE)
```

---

## 5. 统计/分析方法清单

### 5.1 Procrustes 分析

| 参数 | 值 |
|------|-----|
| **函数** | `vegan::procrustes()` |
| **目标矩阵** | 代谢组 PCA 坐标 (PC1, PC2) |
| **拟合矩阵** | 微生物组 PCoA 坐标 (Axis1, Axis2) |
| **对称模式** | `symmetric = TRUE`（双向对称 Procrustes） |
| **输出指标** | m² (sum of squared residuals) |

### 5.2 PROTEST 置换检验

| 参数 | 值 |
|------|-----|
| **函数** | `vegan::protest()` |
| **零假设** | 两矩阵的关联等于随机预期（m² = 1） |
| **置换次数** | 10,000 次 |
| **对称模式** | `symmetric = TRUE` |
| **输出指标** | m², P 值 |
| **随机种子** | `set.seed(123)` |

### 5.3 降维方法汇总

| 数据类型 | 方法 | 距离度量 | 函数 |
|----------|------|----------|------|
| 代谢组 | PCA | 欧氏距离（隐式） | `prcomp()` |
| 微生物组 | PCoA | Bray-Curtis | `vegdist()` + `cmdscale()` |

---

## 6. 可视化编码要点

### 6.1 颜色映射

```r
col_palette <- c(
  "0%"   = "#E64B35FF",  # 红色
  "25%"  = "#4DBBD5FF",  # 青色
  "50%"  = "#00A087FF",  # 绿色
  "75%"  = "#3C3C3CFF",  # 深灰
  "100%" = "#F39B7FFF"   # 橙色
)
```
> 调色板风格：NPG (Nature Publishing Group) 风格

### 6.2 形状映射

| 数据类型 | 形状 | ggplot2 编码 |
|----------|------|--------------|
| Exudates (代谢组) | 实心圆 | `shape = 21` (circle with fill) |
| Microbiota (微生物组) | 实心三角 | `shape = 24` (triangle with fill) |

### 6.3 箭头编码

- **连接**：同一样本的代谢组坐标 → 微生物组（旋转后）坐标
- **箭头长度**：Procrustes 残差（匹配程度的度量）
- **箭头样式**：`arrow(length = unit(0.12, "cm"), type = "closed")`
- **透明度**：`alpha = 0.7`

### 6.4 Panel B（残差箱线图）

- X 轴：处理组（0% ~ 100%）
- Y 轴：箭头长度（Residual Distance）
- 可视化元素：箱线图 + jitter 散点
- 图例：隐藏（与 Panel A 共用颜色编码）

### 6.5 统计标注

图的副标题显示：
```
m² = 0.47, P < 0.001 (PROTEST, 10,000 permutations)
```

### 6.6 输出规格

| 格式 | 尺寸 | 分辨率 |
|------|------|--------|
| PDF | 11 × 5.2 inch | 矢量 |
| TIFF | 11 × 5.2 inch | 600 dpi, LZW 压缩 |

---

## 7. 使用的关键 R 包/函数

### 7.1 统计分析

| 包 | 函数 | 用途 |
|----|------|------|
| **vegan** | `vegdist()` | Bray-Curtis 距离计算 |
| **vegan** | `procrustes()` | Procrustes 旋转与对齐 |
| **vegan** | `protest()` | PROTEST 置换检验 |
| **stats** | `prcomp()` | 主成分分析 (PCA) |
| **stats** | `cmdscale()` | 经典多维标度法 (PCoA) |

### 7.2 可视化

| 包 | 函数 | 用途 |
|----|------|------|
| **ggplot2** | `geom_point()`, `geom_segment()` | 散点与箭头 |
| **ggplot2** | `geom_boxplot()`, `geom_jitter()` | 残差分布 |
| **patchwork** | `wrap_plots()`, `plot_layout()` | 多面板拼接 |
| **grid** | `arrow()` | 箭头样式控制 |

---

## 8. 潜在风险点与注意事项

> [!NOTE]
> 以下为静态代码审查发现的潜在问题，未实际运行验证。

### 8.1 样本量与统计检验力

- 每组仅 **3 个生物学重复**，共 15 个样本
- Procrustes 分析对样本量较敏感，n = 15 时置换检验的统计检验力有限

### 8.2 降维空间的选择

- 仅使用 **前 2 个维度** 进行 Procrustes 分析
- 如果 PC3/Axis3 解释了较大比例的方差，可能丢失重要信息
- **建议**：在正文或补充材料中报告累计方差解释率

### 8.3 距离度量的一致性

- 代谢组使用 **欧氏距离**（PCA 隐式）
- 微生物组使用 **Bray-Curtis 距离**（组成型数据）
- 两种距离度量的生态学含义不同，但 Procrustes 比较的是样本在各自降维空间中的相对位置，而非原始距离值

### 8.4 m² 解释

- `m² = 0.47` 表示两组数据约 **53% 的变异未被匹配**
- ToDO 文档建议在 Discussion 中讨论未解释变异的可能来源

---

## 9. 脚本运行信息

> [!WARNING]
> 当前为静态分析，未实际执行脚本。以下为预期输出。

### 预期输出文件

| 文件 | 路径 |
|------|------|
| `Fig8_combined.pdf` | `./figure8/` |
| `Fig8_combined.tiff` | `./figure8/` |

### 包版本要求

- R ≥ 4.4.3
- vegan, ggplot2, patchwork

---

## 10. 参考文献

1. Peres-Neto PR, Jackson DA (2001) How well do multivariate data sets match? The advantages of a Procrustean superimposition approach over the Mantel test. *Oecologia* 129: 169–178.
2. Oksanen J et al. (2022) vegan: Community Ecology Package. R package version 2.6-4.
