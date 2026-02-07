# Figure 9 分析报告

> **生成日期**: 2026-02-07  
> **脚本**: `figure9_make_fig9.R` (975 行)  
> **R 版本**: 4.4.3

---

## 1. 文件清单

| 类别 | 文件路径 |
|------|----------|
| **主脚本** | `figure9/figure9_make_fig9.R` |
| **ToDO** | `figure9/ToDO_fig9.md` |
| **数据输入** | `data/完整数据-分泌物.csv`（代谢组）, `data/完整数据-微生物.csv`（微生物组） |
| **输出图** | `Fig9_main.pdf/.tiff`, `Fig9A_scores.pdf`, `Fig9B_loadings.pdf`, `Fig9C_chord.pdf` |
| **补充图** | `FigS9_splsda_CV.pdf` |
| **表格** | `splsda_tuning_results.csv`, `FigS9_selected_features.csv`, `FigS9_splsda_loadings_nonzero.csv`, `FigS9_chord_links.csv`, `TableS_CV_tuning.csv`, `TableS_chord_labels.csv` |

---

## 2. 图展示内容

### Panel A: sPLS-DA 得分图 (Scores Plot)
展示不同根系分泌物处理组（0%、25%、50%、75%、100%）在 sPLS-DA 模型空间中的分离情况。通过 95% 置信椭圆可视化各组样本聚类程度，回答"不同浓度梯度的根系分泌物添加是否导致代谢组-微生物组联合特征的系统性差异"。

### Panel B1/B2: 特征载荷图 (Loadings Plot)
- **B1 (Comp 1)**: 显示对第一主成分贡献最大的 Top 10 特征（代谢物/菌群），标注各特征在组间分离中的作用方向与强度。
- **B2 (Comp 2)**: 同上，显示第二主成分的关键特征载荷。

### Panel C: Chord 弦图 (Chord Diagram)
展示筛选后的 15 种代谢物与 12 个微生物分类群（门水平）之间的 Spearman 相关性网络。弦的颜色表示相关方向（红=正相关，蓝=负相关），粗细表示相关强度。外环热图为代谢物的 log₂FC 值（25-100% vs 0%）。

---

## 3. 数据来源与样本结构

| 项目 | 内容 |
|------|------|
| **代谢组** | `完整数据-分泌物.csv`，样本列名格式 `YS{ratio}_F-{rep}` |
| **微生物组** | `完整数据-微生物.csv`，样本列名格式 `RDYS_{group}-{rep}` |
| **分组水平** | 5 组：0%、25%、50%、75%、100% |
| **每组样本数** | 3 个生物学重复 |
| **总样本量** | 15 个样本 |
| **数据层级** | 样本级（单点），非均值 |
| **微生物聚合** | 门水平 (Phylum) |

> [!NOTE]  
> 微生物列名 `RDYS_*` 通过函数 `map_rdys_to_ys()` 映射为 `YS*_F-*`，与代谢组样本名对齐。

---

## 4. 数据预处理流程

### 4.1 代谢组
1. **读取**: `data.table::fread()` / `read.csv()` with UTF-8-BOM
2. **缺失/非有限值处理**: 替换为 0
3. **对数变换**: `log₂(x + 1)`
4. **方差过滤**: 剔除 SD = 0 的特征
5. **矩阵转置**: 输出 samples × metabolites 矩阵

### 4.2 微生物组
1. **门级聚合**: 根据 `taxonomy` 列提取门 (`p__`)，`rowsum()` 合并
2. **相对丰度**: 列归一化 (每样本总和为 1)
3. **低丰度过滤**: 保留平均相对丰度 ≥ 10⁻⁴ 的门
4. **最大特征数限制**: 最多保留 200 个门 (按丰度降序)
5. **对数变换**: `log₂(rel + 10⁻⁶)`
6. **方差过滤**: 剔除 SD = 0 的特征

### 4.3 合并
- 取两矩阵的交集样本 (共 15 个)
- 按列拼接为联合 X 矩阵

---

## 5. 统计/分析方法清单

### 5.1 sPLS-DA (Sparse Partial Least Squares Discriminant Analysis)

| 参数 | 值 |
|------|-----|
| 函数 | `mixOmics::splsda()` |
| 最大组件数 | `ncomp_max = 2` |
| 变量保留网格 | `keepX_grid = c(5, 10, 15, 20, 25, 30)` |
| 标准化 | `scale = TRUE` (内部 Z-score) |
| 随机种子 | `set.seed(123)` |

### 5.2 交叉验证与调参

| 参数 | 值 |
|------|-----|
| 函数 | `mixOmics::tune.splsda()` + `perf()` |
| 验证方式 | M-fold CV |
| 折数 | `folds = 3` |
| 重复次数 | `repeats = 50` |
| 距离度量 | `centroids.dist` |
| 评价指标 | BER (Balanced Error Rate) |

> [!IMPORTANT]  
> CV 折数 = 3 是因为每组仅有 3 个样本，无法使用更大折数。

### 5.3 特征选择

- **Top 代谢物数**: 15
- **Top 门类数**: 12
- **选择依据**: 各成分载荷绝对值排序

### 5.4 相关性分析 (Chord 图)

| 参数 | 值 |
|------|-----|
| 方法 | Spearman 秩相关 |
| 多重校正 | BH-FDR (`p.adjust(..., method = "BH")`) |
| |rho| 阈值 | 0.6 |
| q 值阈值 | 0.05 |
| 高亮阈值 | |rho| ≥ 0.8 时加粗 1.6 倍 |

### 5.5 log₂FC 计算

```
log₂FC = log₂[(mean(25-100%) + ε) / (mean(0%) + ε)]
ε = 10⁻⁵
```

---

## 6. 可视化编码要点

### Panel A (Scores)
- **点**: 按处理组着色
- **椭圆**: 95% 置信椭圆 (ncomp ≥ 2 时启用)
- **坐标轴**: Variates (Comp 1 vs Comp 2)

### Panel B (Loadings)
- **条形方向**: 水平条形图
- **排序**: 按载荷绝对值降序
- **显示数**: 各成分 Top 10
- **正/负**: 条形正负表示贡献方向

### Panel C (Chord)
| 元素 | 编码 |
|------|------|
| **扇区** | 代谢物按 Super_Class 着色；门类统一灰色 |
| **弦颜色** | 红色渐变 = 正相关；蓝色渐变 = 负相关 |
| **弦粗细** | 基础 0.6pt，按 (|rho| - 0.6)/(1 - 0.6) 线性加粗 |
| **高亮** | |rho| ≥ 0.8 再加粗 1.6 倍 |
| **外环热图** | log₂FC (25-100% vs 0%)；红 = 上调，蓝 = 下调 |
| **标签** | 短代码 M1-Mn / T1-Tm，对照表见 `TableS_chord_labels.csv` |

---

## 7. 使用的 R 包/函数

| 包 | 核心函数 | 用途 |
|----|----------|------|
| **mixOmics** | `splsda()`, `tune.splsda()`, `perf()`, `plotIndiv()`, `plotLoadings()` | sPLS-DA 建模与可视化 |
| **circlize** | `chordDiagram()`, `circos.*()`, `colorRamp2()` | Chord 弦图 |
| **data.table** | `fread()` | 高效读取大 CSV |
| **grDevices** | `colorRampPalette()`, `hcl.colors()`, `pdf()`, `tiff()` | 颜色与图形输出 |
| **stats** | `cor.test()`, `p.adjust()` | Spearman 相关 + BH 校正 |

---

## 8. 潜在风险点

| 风险 | 状态 | 说明 |
|------|------|------|
| 组均值替代样本点 | ❌ 无 | 脚本使用样本级数据，非均值 |
| n 值错误 | ❌ 无 | 共 5 组 × 3 重复 = 15 样本，与实际一致 |
| citrate 命名 | ⚠️ 需检查 | 脚本第 567-573 行检测是否存在 citrate/citric acid，如正文提及但数据无，需纠正 |
| 小样本 CV | ⚠️ 注意 | 每组仅 3 样本，3-fold CV 可能不稳定，但 50 次重复可部分缓解 |

---

## 9. 产出文件摘要

| 文件名 | 内容 |
|--------|------|
| `Fig9_main.pdf/.tiff` | 主图 (A + B1/B2 + C 组合) |
| `Fig9A_scores.pdf` | 单独的 sPLS-DA 得分图 |
| `Fig9B_loadings.pdf` | 单独的载荷图 (Comp1 + Comp2) |
| `Fig9C_chord.pdf` | 单独的 Chord 弦图 |
| `FigS9_splsda_CV.pdf` | CV 性能图 (tune + perf 结果) |
| `splsda_tuning_results.csv` | 调参结果 (ncomp, keepX, CV 参数) |
| `TableS_CV_tuning.csv` | 详细 CV 错误率网格 |
| `FigS9_selected_features.csv` | 筛选的 15 代谢物 + 12 门类 |
| `FigS9_splsda_loadings_nonzero.csv` | 所有非零载荷特征 |
| `FigS9_chord_links.csv` | Chord 图使用的边表 (rho, p, q) |
| `TableS_chord_labels.csv` | 短代码与全名对照表 |
