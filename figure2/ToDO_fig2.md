Figure 2（代谢物热图 + Top响应）
来自“总结”的优先级提醒
* Figure 2：主要是排版/可读性与叙述一致性（聚类解释准确）。
原图
* Figure 2A（Top50 heatmap）+ 2B（Top50 log2FC排名）。
修改要求（补齐点）
* 修复：字体/拉伸/一致性（Reviewer#1点名）。
* 文字解释必须与聚类结果一致（Reviewer#2点名“50%更像T3/T4”的问题）。
* 图注补：n=3 trees per treatment、标准化方式（TIC、Z-score、log2(x+1)）写清。
分析方法（R 4.4.3）
* 数据预处理：TIC、log2(x+1)、Z-score；Faith’s PD 需根据原始gz文件分析补充。
* Heatmap：ComplexHeatmap。
* Top50排序：按 max(|log2FC|)（或直接用表里已有 log2FC_Processing_*）。
需要导入的数据
* /mnt/data/完整数据-分泌物.csv（15个YS样本列 + Super_Class）。
* （可选增强筛选/注释）/mnt/data/代谢物定性定量结果表.xlsx（POS/NEG 的 Total score/RSD/QC列）。
产出
* Fig2_main.pdf（矢量）+ Fig2_main.tiff（600 dpi）。
* （可选）Fig2_data_top50.csv（50个代谢物ID/Name/superclass）。
关联补充/新增（Supplementary）
* 新增图 S｜代谢组 PCA loadings/biplot：prcomp + top loadings条形图/双标图；数据：完整数据-分泌物.csv（可结合xlsx过滤高RSD特征）；产出：FigS_metabolome_PCA_loadings.pdf + PCA_loadings_table.csv。