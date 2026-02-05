Figure 6（关键门/菌门丰度轨迹）
来自“总结”的优先级提醒
* Figure 6：必须补 SD/SE + 显著性检验（Reviewer#1 明确点名缺 stdev/statistics）。
原图
* Figure 6（bar plot：relative abundance vs RS）。
修改要求（必须补齐）
* 必须加 均值±SD + 组间检验（每个门一个检验，BH校正）。
* 图注里“site-level means”不严谨：改为“treatment-level mean across trees (n=3)”。
分析方法（R 4.4.3）
* taxonomy解析→phylum聚合→相对丰度。
* 每个 phylum 做 KW + Dunn；或只对展示的代表门做检验（更稳妥）。
* 绘图：point-range（点+误差线）或 bar+errorbar（点更清楚）。
需要导入的数据
* /mnt/data/完整数据-微生物.csv。
产出
* Fig6_main.pdf + Fig6_main.tiff。
* Fig6_phylum_stats.csv（统计表）。