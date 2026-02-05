Figure 7（代表代谢物的非线性响应 + segmented regression）
来自“总结”的优先级提醒
* Figure 7：必须补 segmented regression 方法 + 断点CI；关键确认/修正“Citrate”变量来源。
原图
* Figure 7（log2FC vs RS + 分段回归）。
修改要求（必须补齐）
* Methods/图注写清：segmented regression（Muggeo类）+ 断点 95%CI。
* Reviewer#1点名“where is citrate?”；但数据表未检到 citric acid/citrate（仅有 Triethylcitrate/Citrulline 等）。
* 执行要求：先核对“Citrate”是否真实存在且注释可靠；否则 Figure 7/9 里的 citrate 必须更正（致命一致性问题）。
* 图形：每条曲线建议显示样本点/均值，不要只有线；统一字体。
分析方法（R 4.4.3）
* 用 /完整数据-分泌物.csv 计算 log2FC（相对 Control 均值）。
* 断点回归：segmented；输出 breakpoint ± CI。
* 绘图：ggplot2 + 拟合线 + CI ribbon。
需要导入的数据
* /mnt/data/完整数据-分泌物.csv。
* （可选增强注释/筛选）/mnt/data/代谢物定性定量结果表.xlsx。
* （额外需要）原文 Figure 7 里选取的代谢物名单（从原图/原文提取，或用 Name 匹配自动找）。
产出
* Fig7_main.pdf + Fig7_main.tiff。
* TableS_breakpoints.csv（新增：每个代谢物 breakpoint 与 CI）。
* FigS7_breakpoint_diagnostics.pdf（可选：残差/拟合诊断）。
关联补充/新增（Supplementary）
* 新增4（Table S）：断点回归 breakpoint ± 95%CI（每条代表代谢物/指标一行）。
* 新增表 S｜断点回归 breakpoint ± 95%CI：方法 segmented + confint；数据 完整数据-分泌物.csv（及代表代谢物列表）；产出 TableS_breakpoints.csv。
