# Figure 7 报告

> **生成日期**：2026-02-07  
> **脚本路径**：`figure7/fig7_segmented_regression.R`

---

## 一、相关文件清单

| 类型 | 文件名 | 说明 |
|------|--------|------|
| 主脚本 | `fig7_segmented_regression.R` | Muggeo 分段回归主脚本 |
| 辅助脚本 | `figure7 关键分泌物 log2FC 趋势图（Top50.R` | 早期趋势图脚本（未用于最终图） |
| TODO 文件 | `ToDO_fig7.md` | 修改要求与分析说明 |
| 输入数据 | `data/完整数据-分泌物.csv` | 根系分泌物代谢组原始丰度 |
| 输出图 | `Fig7_main.pdf` / `Fig7_main.tiff` | 分段回归拟合图 |
| 输出表 | `TableS_breakpoints.csv` | 断点估计值 ± 95% CI 汇总 |

---

## 二、本图展示内容

### Panel 全览（facet_wrap，共 6 个代谢物）

Figure 7 展示 **代表性根系分泌物** 对根系切断程度 (RS, %) 的 **非线性响应**，每个子图为一种代谢物：

- **X 轴**：根系切断强度 RS（0%=Control, 25%, 50%, 75%, 100%）
- **Y 轴**：相对于 Control 均值的 log₂Fold Change (log₂FC)
- **拟合曲线**：分段回归（Muggeo）+ 95% CI ribbon
- **红色虚线 + 浅红色带**：断点 (breakpoint) 位置及其 95% CI

**本图回答的问题**：
1. 代表代谢物的释放量随 RS 增加呈 **阈值响应** 还是 **线性响应**？
2. 若存在阈值，**断点位置**（breakpoint）在哪？对应约多大程度的根系损伤？
3. 关键化合物（如 Daidzein、Triethylcitrate）的响应模式是否一致？

---

## 三、数据来源与样本结构

| 项目 | 说明 |
|------|------|
| 数据文件 | `data/完整数据-分泌物.csv` |
| 样本命名规则 | `YS{RS}_F-{Rep}`，如 `YS0_F-1`～`YS100_F-3` |
| 组别（RS 水平） | 0%（Control）、25%、50%、75%、100%，共 5 组 |
| 重复数 | 每组 n = 3 |
| 分析单位 | **样本级**（非 tree-level 均值）；图中同时展示单点（jittered）与组均值 ± SE |
| 所选代谢物 | 固定指定 Triethylcitrate、Daidzein；其余 4 种通过断点筛选自动选出 |

> ⚠️ **Citrate 问题**：原数据中无精确命名为 `Citrate` 的条目，仅有 `Triethylcitrate`、`Citrulline` 等类似物。脚本已将 `Triethylcitrate` 纳入代替。

---

## 四、数据预处理流程

1. **读取 CSV**：`read.csv(..., fileEncoding = "UTF-8-BOM")`，保留原始列名。  
2. **匹配样本列**：正则 `^YS[0-9]+_F-[0-9]+$` 识别样本列并提取 RS 与重复编号。  
3. **计算 Control 均值**：对 RS == 0 的列取行均值作为每个化合物的基线。  
4. **计算 log₂FC**：`log₂(Intensity / ControlMean)`。  
5. **筛选代表代谢物**：
   - 优先纳入固定化合物（Triethylcitrate, Daidzein）；
   - 对其余化合物逐一拟合分段回归，选取断点接近 75% 且 |max log₂FC| > 1 的化合物，最终取 6 条。

---

## 五、统计/分析方法清单

| 方法 | 详细说明 |
|------|----------|
| **分段回归** | `segmented::segmented()` (Muggeo 类)，单断点模型 |
| 初始断点猜测 | `psi = c(25, 50, 75)` |
| 断点 95% CI | `confint(seg_fit, level = 0.95)` |
| 线性回归 fallback | 若分段回归不收敛，使用普通线性回归 `lm(log2FC ~ RS)` |
| 模型评价指标 | R²_adj、AIC、Davies test P 值（`segmented::pscore.test()`） |
| 预测置信带 | `predict(model, se.fit = TRUE)` + t 分布置信区间 |

> **置换次数/多重校正**：本脚本未涉及多次假设检验多重校正（单化合物独立建模）。

---

## 六、可视化编码要点

| 图层 | 说明 |
|------|------|
| **灰色散点**（jitter） | 原始 n=3 样本点，`position_jitter(width=1.8)` |
| **黑色实心点** | 组均值 (`stat_summary(fun=mean)`) |
| **黑色误差棒** | 组均值 ± SE (`stat_summary(fun.data=mean_se)`) |
| **蓝色曲线** | 分段/线性回归拟合线 |
| **蓝色浅带** | 95% 预测置信区间 ribbon |
| **红色虚线** | 断点位置 |
| **红色浅带** | 断点 95% CI 区间 |
| **统计文字标注** | 左上角显示 BP、CI、R²、P 值 |
| 分面 | `facet_wrap(~Compound, scales = "free_y", ncol = 3)` |
| 字体 | `base_family = "sans"`，`strip.text` 加粗 |

---

## 七、使用的关键 R 包/函数

| 包 | 主要函数 | 用途 |
|----|----------|------|
| **segmented** | `segmented()`, `confint()`, `pscore.test()` | 分段回归与断点估计 |
| **ggplot2** | `ggplot()`, `geom_point()`, `geom_ribbon()`, `facet_wrap()`, `ggsave()` | �bind�图与输出 |
| **stats** | `lm()`, `predict()`, `aggregate()` | 线性回归与汇总 |
| **grDevices** / **ragg** | `pdf()`, `agg_tiff()` | 高分辨率输出 |

---

## 八、断点回归结果摘要（TableS_breakpoints.csv）

| Compound | Breakpoint (%) | 95% CI | R²_adj | P 值 |
|----------|----------------|--------|--------|------|
| Triethylcitrate | 25.0 | [−16.9, 25.0] | 0.13 | 0.908 |
| Daidzein | 75.0 | [68.7, 75.0] | 0.88 | **0.033** |
| (4r)-4-hydroxy-5,5-dimethyloxolan-2-one | 75.0 | [55.5, 75.0] | 0.83 | 0.246 |
| N-Acetylalanine | 75.0 | [38.7, 75.0] | 0.17 | 0.813 |
| trans-Cinnamic acid | 75.0 | [41.8, 75.0] | 0.39 | 0.264 |
| vanillic acid | 75.0 | [57.4, 75.0] | 0.46 | 0.316 |

> **解读**：Daidzein 的断点显著（P = 0.033），其余化合物 P > 0.05 但断点估计依然集中于 75%。

---

## 九、潜在风险点

1. **Triethylcitrate 替代 Citrate**：数据中无精确 Citrate 条目，已改用 Triethylcitrate，需在论文中说明或删除。
2. **CI 下限负值**：Triethylcitrate 断点 CI 下限为 −16.9%，物理不可解释；论文呈现时应截断至 [0, 100]。
3. **样本量偏小 (n = 3)**：断点估计 CI 较宽，统计效力有限。
4. **断点全落边界**：多数断点逼近 75% 或 25%，可能受初始猜测影响，需核查模型收敛。

---

## 十、产出文件列表

| 文件 | 说明 |
|------|------|
| `Fig7_main.pdf` | 矢量版本主图 |
| `Fig7_main.tiff` | 600 dpi 位图版本 |
| `TableS_breakpoints.csv` | 断点估计结果表 |
| `FigS7_breakpoint_diagnostics.pdf` | （可选）残差诊断图（需设置 `make_diagnostics <- TRUE`） |
