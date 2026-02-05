Figure 5（网络：代谢物网络 + 微生物网络）
来自“总结”的优先级提醒
* Figure 5：网络图必须补图例/拓扑参数/统计说明，并做阈值敏感性（建议放补充）。
原图
* Figure 5A（exudate network |ρ|>0.7）+ 5B（microbe network Louvain）。
修改要求（必须补齐）
* 网络图可读性（字体、只标注 top hubs）。
* 图注补齐：节点颜色含义、模块划分（Louvain）、边正负、阈值、FDR标准。
* 必须补网络拓扑参数（节点/边/平均度/模块度Q等）。
* 建议新增：阈值敏感性（|ρ|≥0.7 vs 0.5）对照。
分析方法（R 4.4.3）
* 相关：Spearman + BH-FDR（Hmisc::rcorr / psych::corr.test）。
* 构网：igraph；模块：cluster_louvain。
* 绘图：ggraph。
* 敏感性：两套阈值下输出 hub overlap / 模块度对比。
需要导入的数据
* /mnt/data/完整数据-分泌物.csv（代谢矩阵）。
* /mnt/data/完整数据-微生物.csv（建议先聚合到 phylum/genus 并筛 Top N）。
产出
* Fig5_main.pdf + Fig5_main.tiff。
* Fig5_network_topology.csv（拓扑指标表）。
* FigS5_threshold_sensitivity.pdf（新增：0.7 vs 0.5 对照图）。
* FigS5_threshold_sensitivity_table.csv（新增：对照数值）。
关联补充/新增（Supplementary）
* 新增5（FigS + Table S）：网络 |ρ|≥0.7 vs 0.5 的敏感性对照 + hub一致性/模块度等。
* 新增图/表 S｜网络阈值敏感性：方法 两套阈值网络 + hub overlap + 模块度Q/平均度对比；数据 分泌物+微生物（聚合后）；产出 FigS_network_sensitivity.pdf + network_sensitivity_metrics.csv。
* 新增6（FigS）：正负相关边比例（2:1）置换检验零分布图。
* 新增图 S｜正负边比例置换检验：方法 置换产生零分布 + 观测值竖线；数据 来自网络边表（T10输出）；产出 FigS_posneg_permutation.pdf + posneg_test_results.csv。
