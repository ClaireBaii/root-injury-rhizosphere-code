Figure 3（α多样性 + PCoA）
来自“总结”的优先级提醒
* Figure 3：必须补误差线 + 显著性检验细节；建议新增 distance-to-control 量化面板。
原图
* Figure 3A（Chao1/FaithPD/Shannon/Simpson）+ 3B（PCoA）。
修改要求（必须补齐）
* α多样性：加每个样本点 + 均值±SD，并补 Kruskal-Wallis + Dunn(BH) 显著性标注。
* 数字叙述与坐标轴一致（Reviewer#2点名 Chao1 drop 描述可能不对）。
* PCoA：图注写清 PERMANOVA 的 R²、P、置换次数；“100%回撤”建议用 distance-to-control 量化支撑。
分析方法（R 4.4.3）
* 从 counts 计算：Chao1/Shannon/Simpson（phyloseq 或 vegan）。
* PCoA：Bray–Curtis（vegan::vegdist + cmdscale / ape::pcoa）。
* PERMANOVA：vegan::adonis2（permutations=999 或 9999）。
* 显著性：rstatix::kruskal_test + dunn_test + ggpubr::stat_pvalue_manual。
需要导入的数据
* /mnt/data/完整数据-微生物.csv（OTU counts + taxonomy）以及原始测序数据。
产出
* Fig3_main.pdf + Fig3_main.tiff。
* Fig3_alpha_stats.csv（组均值±SD + KW p + Dunn q）。
* （可选）Fig3_pcoa_scores.csv（PCoA坐标+组别）。
关联补充/新增（新panel或Supplementary）
* 新增1（Fig3新panel或FigS）：各处理组到 Control 的 Bray–Curtis distance-to-control（均值±SD+检验）。
* 新增图 S（建议命名：FigS3F 或 FigS10）｜distance-to-control：目的量化“100% RS部分回撤”（C049）；方法 Bray–Curtis distance-to-control（两种定义任选其一，并在补充给另一种一致性）；数据 完整数据-微生物.csv；产出 FigS_distance_to_control.pdf + .tiff + distance_to_control_stats.csv。
* 新增2（FigS）：微生物 rarefied vs non-rarefied 的 PCoA/PERMANOVA 稳健性对照。
* 新增图 S｜rarefied vs non-rarefied：方法 rarefy_even_depth + PCoA + PERMANOVA（两套并排）；数据 完整数据-微生物.csv；产出 FigS_rarefaction_robustness.pdf + permanova_robustness.csv。