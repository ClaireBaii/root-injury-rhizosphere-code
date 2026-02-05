Figure 8（Procrustes：代谢组PCA vs 微生物PCoA）
原图
* Figure 8（m²=0.47, P=0.01）。
修改要求（补齐点）
* 图本身通常不缺统计；但建议在 Discussion 补充未解释部分（~53%）可能来源。
* Phase2 加强：可在补充给 Procrustes 残差/箭头长度分布。
分析方法（R 4.4.3）
* vegan::procrustes + protest(permutations=10000)。
* 代谢组 ordination：PCA（prcomp）。
* 微生物 ordination：PCoA（Bray–Curtis）。
需要导入的数据
* /mnt/data/完整数据-分泌物.csv；/mnt/data/完整数据-微生物.csv。
产出
* Fig8_main.pdf + Fig8_main.tiff。
* （可选）FigS8_procrustes_residuals.pdf。