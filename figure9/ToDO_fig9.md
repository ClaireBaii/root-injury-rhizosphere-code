Figure 9（sPLS-DA + loadings + chord）
来自“总结”的优先级提醒
* Figure 9：必须补 sPLS-DA CV细节 + 弦图可读性 + 关键配对强调。
原图
* Figure 9A（sPLS-DA score）+ 9B（loadings）+ 9C（chord）。
修改要求（必须补齐）
* sPLS-DA：补 ncomp、keepX、CV folds/repeats、变量数如何定。
* Chord：字体放大、图例写清；正文/图注强调关键配对。
* 同 Figure 7：若正文写了 citrate 但表里不存在，需要纠正。
分析方法（R 4.4.3）
* mixOmics::splsda + tune.splsda + perf（输出CV性能）。
* chord：筛选后的 15×12 相关矩阵；circlize 绘制。
* 相关检验：Spearman + BH-FDR（|ρ|>0.6、q<0.05按稿子）。
需要导入的数据
* /mnt/data/完整数据-分泌物.csv（代谢矩阵）。
* /mnt/data/完整数据-微生物.csv（建议聚合到 phylum 或 genus，与原文一致）。
产出
* Fig9_main.pdf + Fig9_main.tiff。
* FigS9_splsda_CV.pdf（新增：CV性能/误差率图）。
* FigS9_selected_features.csv（新增：15代谢物+12类群清单）。
* FigS9_chord_links.csv（新增：弦图边表）。
关联补充/新增（Supplementary）
* 新增7（FigS）：sPLS-DA perf/tuning 结果图（BER/错误率 vs keepX/ncomp）。
* 新增图 S｜sPLS-DA CV性能图：对应审稿点 C055；方法 tune.splsda + perf；数据 分泌物矩阵 + 微生物聚合矩阵；产出 FigS_splsda_CV.pdf + splsda_tuning_results.csv。
