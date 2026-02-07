# Figure 2 技术报告

> **生成时间**: 2026-02-07  
> **目标**: 为论文 Figure 2 提供完整的方法学描述与可复现性文档

---

## 1. 文件清单

| 类型 | 文件名 | 路径 |
|------|--------|------|
| **主脚本** | `fig2_make_main.R` | `./figure2/` |
| **辅助脚本** | `fig2_helpers.R` | `./figure2/` |
| **补充图脚本** | `figS_metabolome_PCA_loadings.R` | `./figure2/` |
| **任务说明** | `ToDO_fig2.md` | `./figure2/` |
| **数据输入** | `完整数据-分泌物.csv` | `./data/` |
| **可选数据** | `代谢物定性定量结果表.xlsx` | `./data/` |
| **运行日志** | `run_fig2.log` | `./figure2/` |
| **PDF 输出** | `Fig2_main.pdf` | `./figure2/` |
| **TIFF 输出** | `Fig2_main.tiff` (600 dpi) | `./figure2/` |
| **中间数据** | `Fig2_data_top50.csv` | `./figure2/` |

---

## 2. 图形内容解读

### Panel A — 热图 (Heatmap)

**展示内容**: Top 50 差异代谢物的丰度模式（TIC 归一化 → log2(x+1) → Z-score 行标准化后）

**回答的科学问题**: 不同盐胁迫处理 (YS0, YS25, YS50, YS75, YS100) 下，杨树根系分泌物中哪些代谢物表现出最显著的响应？各处理组之间的代谢谱模式如何？

### Panel B — 右侧注释条形图 (Max log2FC Bar)

**展示内容**: Top 50 代谢物的最大 log2 Fold Change（相对于对照组 YS0），按绝对值降序排列

**回答的科学问题**: 哪些代谢物的响应程度最大？它们是上调还是下调？在哪个处理组的响应最强？

---

## 3. 数据来源与样本结构

| 项目 | 描述 |
|------|------|
| **数据文件** | `./data/完整数据-分泌物.csv` |
| **样本列格式** | `YS{处理}_F-{重复号}` (如 `YS0_F-1`, `YS100_F-3`) |
| **处理组** | YS0, YS25, YS50, YS75, YS100 |
| **每组样本数** | **n = 3** (tree-level, 3 棵树/处理) |
| **代谢物注释** | 含 `Name`, `Super_Class` 等列 |
| **绘图单位** | **组均值** (3 个样本取平均后绘制) |

> [!WARNING]
> **风险点**: 热图绘制使用 **组均值 (n=3 → 1)**，而非逐样本(n=15)展示。这意味着图中只有 5 列（5 个处理组），而不是 15 列（15 个样本）。
> 
> 请确保图注中明确说明 "values represent group means (n = 3 trees per treatment)"。

---

## 4. 数据预处理流程

```mermaid
flowchart LR
    A[原始峰面积] --> B[TIC 归一化]
    B --> C["log2(x+1) 转换"]
    C --> D["Z-score 行标准化"]
    D --> E[Top 50 筛选]
    E --> F[热图绘制]
```

### 详细步骤

1. **TIC 归一化 (Total Ion Current normalization)**
   - 公式: `x_norm = (x / TIC) × 10^6`
   - TIC = 每个样本所有代谢物峰面积之和
   - 目的: 校正不同样本间总信号量差异

2. **log2(x+1) 转换**
   - 公式: `x_log = log2(x_norm + 1)`
   - "+1" 是为了避免零值导致的 log(0) 无穷大问题
   - 目的: 压缩数据动态范围，使分布更接近正态

3. **组均值计算**
   - 对每个处理组内的 3 个样本取算术平均
   - 结果: 从 15 列（样本）变为 5 列（处理组）

4. **Z-score 行标准化**
   - 公式: `z = (x - μ_row) / σ_row`
   - 按代谢物（行）进行，使每行均值=0、标准差=1
   - 目的: 便于跨代谢物比较相对高低模式

5. **Top 50 筛选**
   - 计算每个代谢物在各非对照组 vs YS0 的 log2FC
   - 取 `max(|log2FC|)` 作为排序依据
   - 选取排名前 50 的代谢物绘图

6. **缺失值处理**
   - `NA` 值在原始矩阵中被替换为 `0`
   - Z-score 计算时排除标准差为 0 的行（常量行）

---

## 5. 统计/分析方法清单

| 方法 | 用途 | 关键参数 |
|------|------|----------|
| **TIC 归一化** | 样本间信号校正 | scale_factor = 10^6 |
| **log2 Fold Change** | 差异倍数计算 | 公式: log2((Treatment+1)/(YS0+1))，基于 TIC 后均值 |
| **Z-score** | 行标准化热图展示 | center = TRUE, scale = TRUE (对行操作) |

> [!NOTE]
> 本图 **无** 显式统计检验（如 ANOVA、Kruskal-Wallis）。Top 50 的筛选标准是 max(|log2FC|) 的描述性排序，而非基于 p 值或 FDR。

---

## 6. 可视化编码要点

### 颜色映射

| 元素 | 颜色方案 | 说明 |
|------|----------|------|
| **热图** | `#2166AC` (蓝) — `white` — `#B2182B` (红) | Z-score 范围 [-2, 0, 2]，使用 `circlize::colorRamp2` |
| **Super_Class** | Okabe-Ito 色盲友好色板 | 9 色循环: `#E69F00`, `#56B4E9`, `#009E73`, ... |
| **Treatment** | 5 色方案 | `#4053d3` (蓝), `#ddb310` (金), `#b51d14` (红), `#00beff` (青), `#fb49b0` (粉) |
| **log2FC 条形图** | 红/蓝 | 上调 (≥0) 使用 `#B2182B`；下调 (<0) 使用 `#2166AC` |

### 布局与注释

- **左侧注释**: Super_Class（代谢物分类）
- **右侧注释**: Max_Treatment（响应最强的处理组）+ Max log2FC 条形图
- **顶部注释**: Treatment（处理组，顺序 YS0→YS100）
- **聚类**: 行不聚类 (`cluster_rows = FALSE`)，按 max(|log2FC|) 降序排列；列聚类 (`cluster_columns = TRUE`)
- **字体**: 行名 8pt，列名 9pt
- **图尺寸**: PDF/TIFF 均为 15×9 英寸

---

## 7. 使用的关键 R 包/函数

### 统计与数据处理

| 包/函数 | 用途 |
|---------|------|
| `base::rowMeans()` | 计算组均值 |
| `base::scale()` | Z-score 标准化 |
| `base::log2()` | 对数变换 |
| `stats::sweep()` | TIC 归一化 |

### 可视化

| 包 | 用途 |
|----|------|
| **ComplexHeatmap** | 核心热图绑定，支持复杂注释 |
| **circlize** | 颜色渐变函数 `colorRamp2()` |
| **grid** | 底层图形排版、字体控制 |

---

## 8. 运行信息 (来自 run_fig2.log)

```
Figure 2 预处理：TIC → log2(x+1) → Z-score（行）
每组样本数：YS0=3, YS100=3, YS25=3, YS50=3, YS75=3（按 Treatment 取均值绘图）
输出：.../figure2/Fig2_main.pdf
输出：.../figure2/Fig2_main.tiff
```

### 产生的文件

| 文件 | 说明 |
|------|------|
| `Fig2_main.pdf` | 矢量图，可编辑 |
| `Fig2_main.tiff` | 600 dpi 位图，投稿用 |
| `Fig2_data_top50.csv` | Top 50 代谢物数据表（含 max_abs_log2FC, max_Treatment 等） |

---

## 9. 潜在风险与注意事项

> [!CAUTION]
> **关键风险点**
> 
> 1. **组均值替代样本点**: 热图显示的是 5 个组均值（n=3→1），而非 15 个独立样本。这可能隐藏组内变异信息。
> 
> 2. **Top 50 无统计检验**: 当前筛选标准仅基于 log2FC 大小，未进行多组比较的显著性检验（如 ANOVA + Tukey HSD）或多重校正。
> 
> 3. **缺失值处理**: NA 被替换为 0，可能对低丰度代谢物的 log2FC 造成偏倚（分母加 1 可缓解但不完全消除）。

---

## 10. 补充图说明 (Fig S — PCA Loadings)

脚本 `figS_metabolome_PCA_loadings.R` 生成 PCA 双标图与 Top loadings 条形图：

- **预处理**: TIC → log2(x+1)（无 Z-score）
- **PCA 方法**: `prcomp(center = TRUE, scale. = TRUE)`
- **产出**: `FigS_metabolome_PCA_loadings.pdf`, `PCA_loadings_table.csv`
