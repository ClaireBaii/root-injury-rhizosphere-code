Figure 4（门水平热图）
原图
* Figure 4（phylum heatmap，提“两个模块”）。
修改要求（补齐点）
* 图例/字体放大；图注列出“关键门”（Reviewer#2点名：别只写“~27 phyla”）。
* 若主文放不下：可把“完整版热图”放补充，主文保留简化版（Reviewer#3建议）。
分析方法（R 4.4.3）
* taxonomy 解析到 phylum → 聚合相对丰度。
* 热图：ComplexHeatmap（统一字体/注释模块）。
需要导入的数据
* /mnt/data/完整数据-微生物.csv（taxonomy字段需解析）。
产出
* Fig4_main.pdf + Fig4_main.tiff。
* （可选补充）FigS4_full_phylum_heatmap.pdf。
