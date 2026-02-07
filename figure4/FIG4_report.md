# Figure 4 结构化报告：门水平热图

---

## 1. 相关文件列表

| 文件类型 | 路径 |
|---------|------|
| 主脚本 | `figure4/figure4 单样本分组热图.R` |
| TODO 文档 | `figure4/ToDO_fig4.md` |
| 数据输入 | `../data/完整数据-微生物.csv` |
| 输出 (主图) | `Fig4_main.pdf`, `Fig4_main.tiff` |
| 输出 (补充) | `FigS4_full_phylum_heatmap.pdf` |
| 辅助表格 | `Fig4_key_phyla.csv`, `Fig4_key_phyla.txt` |

---

## 2. 图的内容与目的

### Panel A（主图 Fig4_main）

展示 **27 个关键门（Key Phyla）** 在 **5 组处理间的丰度模式差异**。  
- 横轴：5 个处理组（RS0–RS100），代表不同根际土比例。  
- 纵轴：按平均相对丰度筛选的 Top 27 门。  
- 行被划分为 **Module 1 / Module 2**，反映聚类后的丰度响应模块。

> **科学问题**：不同根际土比例如何改变关键微生物门的群落组成？

### Panel B（补充图 FigS4）

展示 **所有非零门 × 所有样本** 的完整热图，供审稿人查看细节。

---

## 3. 数据来源与样本结构

| 项目 | 描述 |
|-----|------|
| 数据来源 | `../data/完整数据-微生物.csv` |
| 样本命名规则 | `RDYS_1-1`, `RDYS_1-2`, ... `RDYS_5-3` （处理×重复） |
| 处理组 | 5 组：RS0, RS25, RS50, RS75, RS100 |
| 每组重复 | n = 3（推断自样本命名 `-1`, `-2`, `-3`） |
| 主图单位 | 处理水平均值（n = 3 平均） |
| 补充图单位 | 单个样本（共 15 列） |

> **注意**：样本分组与重复数从样本命名解析，脚本假设 `RDYS_X-Y` 格式；若原始数据结构不同可能需核实。

---

## 4. 数据预处理流程

1. **读取原始数据**  
   使用 `data.table::fread()` 或 `read.csv()` 读取 CSV。

2. **Taxonomy 解析**  
   正则提取 `p__PhylumName`，未匹配归为 "Unclassified"。

3. **聚合至门水平**  
   `rowsum()` 按 Phylum 汇总 counts。

4. **相对丰度计算**  
   ```r
   rel_mat <- sweep(counts_mat, 2, colSums(counts_mat), "/")
   ```

5. **筛选 Top N 门**  
   取平均相对丰度最高的 27 个门（参数 `top_n_main`）。

6. **处理水平均值**  
   `rowMeans()` 计算每个处理的 n=3 样本均值。

7. **行方向 Z-score 标准化**  
   ```r
   mat_main <- t(scale(t(mat_rel)))
   ```
   对每门进行行标准化（均值=0，SD=1），突出跨处理的变化趋势。

8. **层次聚类 & 模块划分**  
   - 距离：欧氏距离（`dist(z)`）
   - 聚类方法：Ward.D2
   - 分割：`cutree(hc, k = 2)` → Module 1 / Module 2

---

## 5. 统计/分析方法清单

| 方法 | 细节 |
|------|------|
| 层次聚类 | 行方向 Z-score 矩阵 → 欧氏距离 → Ward.D2 → k=2 |
| 模块划分 | `cutree()` 分为 2 个模块 |
| 无假设检验 | 本图不涉及 ANOVA/Kruskal-Wallis 等统计检验 |

> **注意**：热图仅做可视化与聚类，未对门丰度进行差异检验。

---

## 6. 可视化编码要点

| 元素 | 编码规则 |
|------|---------|
| **颜色映射（热图）** | 蓝-白-红渐变（`colorRamp2(c(-1.5, 0, 1.5))`）|
| **列注释（Treatment）** | 5 色：RS0=#1B9E77, RS25=#D95F02, RS50=#7570B3, RS75=#E7298A, RS100=#66A61E |
| **行注释（Module）** | Module 1=#E41A1C（红）, Module 2=#377EB8（蓝） |
| **行分割** | 按模块分割，模块间隔 2 mm |
| **字体大小** | 行名 12 pt, 列名 13 pt, 图例 11-12 pt |
| **输出规格** | PDF 7×9 in；TIFF 600 dpi LZW 压缩 |

---

## 7. 关键 R 包与函数

| 包 | 功能 |
|---|------|
| `ComplexHeatmap` | 热图绑定、注释、分割 |
| `circlize` | `colorRamp2()` 颜色映射 |
| `grid` | 图形参数 (`gpar`) |
| `data.table` | 高效读取 CSV（可选） |

关键函数：
- `Heatmap()` – 创建热图对象
- `HeatmapAnnotation()` / `rowAnnotation()` – 添加列/行注释
- `hclust()` / `cutree()` – 聚类与模块划分
- `scale()` – Z-score 标准化

---

## 8. 风险与建议

| 风险点 | 说明 |
|--------|-----|
| **主图使用组均值** | 主图列为处理均值（n=3），非单样本点；补充图展示全部 15 个样本 |
| **无统计检验** | 热图本身不包含差异显著性分析；如需加星号需额外做 ANOVA/Wilcoxon |
| **门名称解析依赖格式** | 若 taxonomy 字段格式变化（无 `p__`），可能误分类 |

---

## 9. 产出文件汇总

| 文件 | 说明 |
|------|------|
| `Fig4_main.pdf` | 主图（7×9 in，矢量） |
| `Fig4_main.tiff` | 主图（600 dpi，出版） |
| `FigS4_full_phylum_heatmap.pdf` | 补充图（11×10 in） |
| `Fig4_key_phyla.csv` | 关键门列表 + 平均相对丰度 |
| `Fig4_key_phyla.txt` | 关键门名称纯文本 |
