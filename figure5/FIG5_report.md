# Figure 5 结构化报告

> **文件:** `FIG5_report.md`  
> **生成日期:** 2026-02-07

---

## 1. 文件清单

| 类型 | 路径/文件名 |
|------|-------------|
| 主脚本 | `fig5_make_all.R` (519 行) |
| TODO | `ToDO_fig5.md` |
| 数据输入 | `Fig5A_exudate_edges_thr0.8.csv`, `Fig5A_exudate_nodes_thr0.8.csv` (分泌物); `Fig5B_microbe_edges_thr0.6.csv`, `Fig5B_microbe_nodes_thr0.6.csv` (微生物); `Fig5_network_topology.csv` (拓扑表) |
| 辅助脚本 | `figure5 A 分泌物共变网络图.R`, `figure5 B 微生物共现网络图优化脚本.R` (均调用 `fig5_make_all.R`) |
| 输出图 | `Fig5_main.pdf/png`, `Fig5A_module_summary.pdf/png`, `Fig5B_microbiome_network.pdf/png`, `Fig5C_network_metrics.pdf/png` |
| 补充图 | `FigS5A_exudate_network.pdf`, `FigS5B_microbe_network.pdf`, `FigS_posneg_permutation.pdf` |
| 运行日志 | 无独立 log 文件 (无 `run_fig5.log`) |

---

## 2. Panel 分解

### Panel A：分泌物共现网络 (模块汇总)

**展示内容：**  
通过 Louvain 聚类将 1686 个代谢物特征划分为多个模块，取节点数最多的 Top 8 模块（其余归为 "Other"），并显示模块间连边数最多的 Top 15 组模块对关系。节点大小编码模块所含特征数，边宽编码模块间边数量。

**数据来源：**  
- Spearman 相关矩阵 → 阈值筛选 |ρ| ≥ 0.8，BH-FDR < 0.05  
- 来自 `Fig5A_exudate_nodes_thr0.8.csv` / `Fig5A_exudate_edges_thr0.8.csv`  
- 样本层面：全部处理组合并计算；不区分 tree-level；见 topology.csv 中 feature_filter = `max_abs_log2fc_top60`

**分析方法：**  
1. **相关计算:** Spearman's ρ，p < 0.05 后 BH-FDR 校正  
2. **网络构建:** `igraph::graph_from_data_frame`；边权重 = |ρ|  
3. **模块检测:** `igraph::cluster_louvain`  
4. **聚合逻辑:** 统计模块间边数量 (`group_by(from, to) + summarise(edge_count)`)；保留 Top 15  
5. **关于 n 值:** 节点数 = 1686；边数 = 70645；无样本级 n 直接显示

---

### Panel B：微生物共现网络 (稀疏化)

**展示内容：**  
仅保留度数最高的 Top 40 个菌属 (genus)，每个节点仅保留权重最大的 k=2 条边 (sparsify)，标出度数最高的 12 个 hub；节点颜色映射 Louvain 模块（≤6 个模块 + Other）；边颜色区分正/负相关。

**数据来源：**  
- 阈值 |ρ| ≥ 0.6，BH-FDR < 0.05  
- 来自 `Fig5B_microbe_nodes_thr0.6.csv` / `Fig5B_microbe_edges_thr0.6.csv`  
- 数据层面：先聚合到 genus 层级并筛选 Top 40 (prevalence ≥ 0.15)

**分析方法：**  
1. 相关 + FDR 同 Panel A  
2. **子图提取:** 保留度数 Top 40 节点 (`induced_subgraph`)  
3. **边稀疏化:** 每节点取权重 Top-k=2 边 (`incident_edges` + order)  
4. **模块检测:** 在稀疏子图上重新运行 `cluster_louvain`  
5. **hub 标注:** 度数 Top 12 节点显示标签 (`ggrepel::geom_text_repel`)

---

### Panel C：网络拓扑指标表

**展示内容：**  
两行表格，分别对应 Exudate 和 Microbe 网络的拓扑指标：Nodes, Edges, Average Degree, Modularity Q, Positive/Negative Ratio。

**数据来源：**  
- 直接读取 `Fig5_network_topology.csv`

**关键数值 (来自 csv)：**

| Dataset | Nodes | Edges | Avg Degree | Modularity Q | Pos:Neg |
|---------|-------|-------|------------|--------------|---------|
| exudate (thr 0.8) | 1686 | 70645 | 83.80 | 0.598 | 1.59 |
| microbe_genus (thr 0.6) | 40 | 126 | 6.30 | 0.599 | 2.15 |

---

## 3. 数据预处理流程

1. **代谢物矩阵:**  
   - 特征筛选: `max_abs_log2fc_top60` (脚本中可见 feature_filter 标记)  
   - 数值变换: **不确定** — 脚本中未显式提及 log2(x+1) / TIC / Z-score；输入可能为原始或预处理后的峰面积  
2. **微生物矩阵:**  
   - 聚合至 genus → Top 40 (prevalence ≥ 0.15)  
   - 相对丰度转换: **不确定** — 脚本未见 TSS()  
3. **缺失值:** 脚本中均未显式处理 (假设输入已完整)

> [!WARNING]
> 脚本中 **未显式记录** 代谢物及微生物原始矩阵的归一化/标准化步骤。建议在方法中补充此信息。

---

## 4. 统计与分析方法清单

| 方法 | 细节 |
|------|------|
| **相关系数** | Spearman's ρ |
| **多重校正** | Benjamini-Hochberg (FDR < 0.05) |
| **阈值** | Exudate |ρ| ≥ 0.8；Microbe |ρ| ≥ 0.6 |
| **网络构建** | `igraph::graph_from_data_frame`，无向加权 |
| **模块检测** | `igraph::cluster_louvain` (基于模块度优化) |
| **Top-k 稀疏化 (Panel B)** | 每节点保留权重最大的 k=2 条边 |
| **布局算法** | Panel A: circle；Panel B: Fruchterman-Reingold (`layout = "fr"`) |
| **缩放因子** | 布局坐标 ×1.6 (A) / ×1.8 (B) |

---

## 5. 可视化编码要点

| 元素 | 编码 |
|------|------|
| **节点大小** | Panel A = 模块内节点数；Panel B = 稀疏图度数 |
| **节点颜色** | Panel A = 模块 ID (Set2 调色板)；Panel B = 稀疏图 Louvain 模块 (Paired) |
| **边宽** | Panel A = 模块间边数 (range 0.5–3) |
| **边颜色** | Panel B: 红色 = 正相关，蓝色 = 负相关 |
| **标签** | Panel B 使用 `ggrepel`，仅标出度数 Top 12 hub |
| **误差线** | 本图 **无误差线** |
| **显著性** | 通过阈值 + FDR 预先筛选，不在图中标 */* 等 |

---

## 6. R 包与关键函数

| 包 | 主要函数 |
|----|----------|
| `igraph` | `graph_from_data_frame`, `cluster_louvain`, `degree`, `induced_subgraph`, `subgraph.edges`, `incident_edges` |
| `ggraph` | `ggraph`, `create_layout`, `geom_edge_link`, `geom_node_point`, `geom_node_label` |
| `ggplot2` | 底层绑定 |
| `patchwork` | `+`, `/`, `wrap_elements` 多面板拼接 |
| `gridExtra` | `tableGrob`, `arrangeGrob` |
| `ggrepel` | `geom_text_repel` |
| `dplyr` / `tidyr` | 数据整理 |
| `RColorBrewer` | `brewer.pal("Set2", "Paired")` |

---

## 7. 潜在风险点与待确认事项

1. **原始数据标准化未写入脚本：** 无法确认输入矩阵是否已 log2、Z-score 或 TIC 归一化。  
2. **样本层 n 不直接可见：** 网络是基于特征间相关，而非样本点。若需报告样本量，需查阅上游数据文件。  
3. **阈值敏感性：** `ToDO_fig5.md` 提及需 0.7 vs 0.5 对照，但主脚本 config 仅使用 0.8 / 0.6 阈值。若需敏感性分析，需额外调整参数或运行多次。  
4. **正负边比例置换检验：** 脚本中未见执行，ToDO 中标为"新增"。

---

## 8. 输出文件汇总

| 文件 | 说明 |
|------|------|
| `Fig5_main.pdf/png/vector.pdf` | 主图 (220×140 mm, 600 dpi) |
| `Fig5A_module_summary.*` | Panel A 单独 (100×80 mm) |
| `Fig5B_microbiome_network.*` | Panel B 单独 (140×80 mm) |
| `Fig5C_network_metrics.*` | Panel C 单独 (120×45 mm) |
| `FigS5A_exudate_network.pdf` | 完整分泌物网络 (无标签) |
| `FigS5B_microbe_network.pdf` | 完整微生物网络 (Top 10% hub 标签) |
| `Fig5_network_topology.csv` | 拓扑数值表 |
| `FigS5_threshold_sensitivity_table.csv` | 阈值敏感性数值 |

---

*报告完成*
