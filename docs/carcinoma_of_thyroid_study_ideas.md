# carcinoma of thyroid

## Q1
是否有甲状腺癌方面的大的组学数据库, 比如TCGA有没有? 
* Answer
是的，甲状腺癌（thyroid cancer）相关的大规模组学数据库是非常成熟的，而且你提到的 TCGA 就是核心资源之一。下面我帮你系统梳理一下（按“是否主流 + 数据类型 +适用场景”来分层），方便你直接用于科研设计。

⸻

一、核心：TCGA（绝对主力）

✔ The Cancer Genome Atlas (TCGA) — THCA项目

👉 甲状腺癌在 TCGA 中对应项目：TCGA-THCA

数据规模
	•	~500例甲状腺癌患者  ￼
	•	主要是 乳头状甲状腺癌（PTC）为主

多组学类型（非常完整）
	•	WES / mutation
	•	RNA-seq（bulk）
	•	miRNA-seq
	•	DNA methylation（450K）
	•	Copy Number Variation (CNV)
	•	临床数据（分期、生存、突变）

👉 还有：
	•	影像数据（CT/PET）在 TCIA 中配套  ￼
	•	可通过 UCSC Xena / GDC / cBioPortal访问

特点（科研角度）
	•	多组学整合（multi-omics）标准数据源
	•	非常适合：
	•	DEGs分析
	•	生存分析
	•	突变驱动机制（BRAF / RAS）
	•	immune infiltration（CIBERSORT等）

⸻

二、TCGA扩展生态（强烈推荐一起用）

1. Genomic Data Commons (GDC)
	•	TCGA官方数据入口
	•	可下载：
	•	BAM / FASTQ（需授权）
	•	processed counts（公开）

👉 你做pipeline或方法开发，建议直接用这里

⸻

2. UCSC Xena
	•	已整理好的表达矩阵 + clinical
	•	一键下载
	•	适合：
	•	快速建模
	•	machine learning

⸻

3. cBioPortal
	•	强调突变 + pathway
	•	可做：
	•	Oncoprint
	•	co-mutation
	•	survival stratification

⸻

三、GEO（补充数据，尤其单细胞）

✔ Gene Expression Omnibus (GEO)

数据类型
	•	bulk RNA-seq / microarray
	•	scRNA-seq（非常关键）

例如：
	•	scRNA dataset（如 GSE184362）
	•	包含：
	•	原发肿瘤
	•	淋巴转移
	•	正常组织  ￼

👉 用途：
	•	细胞类型解析（T cell exhaustion / CAF）
	•	验证TCGA bulk结果

⸻

四、蛋白组 / 蛋白表达

1. Human Protein Atlas
	•	TCGA转录组 + IHC整合
	•	提供：
	•	prognostic genes
	•	protein expression

👉 很适合：
	•	biomarker筛选
	•	validation

⸻

2. Proteomics（PRIDE / OmicsDI）
	•	DIA / TMT蛋白组
	•	覆盖：
	•	PTC / PDTC / ATC progression  ￼

👉 用途：
	•	multi-omics整合（转录 vs 蛋白）

⸻

五、其他重要数据库（补充）

1. International Cancer Genome Consortium (ICGC)
	•	TCGA的国际版本
	•	甲状腺数据较少，但有补充

⸻

2. SEER
	•	不属于组学
	•	但有：
	•	超大样本（流行病学）
	•	生存分析

⸻

3. 单细胞专用数据库
	•	Human Cell Atlas
	•	TISCH（tumor immune scRNA）

👉 用于：
	•	tumor microenvironment
	•	immune cell解析

⸻

六、从“科研设计”角度怎么组合（重点）

如果你要做论文级 multi-omics thyroid cancer研究，标准组合是：

✅ 推荐组合（经典范式）
	1.	TCGA-THCA
	•	主分析（DEG / survival / mutation）
	2.	GEO（scRNA）
	•	细胞来源验证
	3.	Human Protein Atlas
	•	蛋白层验证
	4.	（可选）Proteomics dataset
	•	multi-omics integration

⸻

七、关键总结（非常重要）

✔ TCGA确实有甲状腺癌，而且是最核心数据源
✔ 项目名：TCGA-THCA
✔ 样本量：~500例（足够做大部分分析）
✔ 数据类型：真正意义上的 multi-omics
✔ 学术界使用率：极高（几乎所有thyroid omics论文都会用）

⸻

如果你下一步要做研究（我可以帮你设计）

你现在的背景（RNA-seq / multi-omics / pipeline开发）其实很适合做：

可以考虑的方向：
	•	immune microenvironment（结合CIBERSORT / scRNA）
	•	ferroptosis / metabolism（你之前在做）
	•	chitin metabolism（你当前研究也能延伸）
	•	multi-omics integration（TCGA + proteomics）

⸻

如果你愿意，我可以帮你直接设计一个**“甲状腺癌 multi-omics 发表级分析框架（含figure layout +方法 pipeline）”**，可以直接冲 Bioinformatics / Cancers / iMeta 这种期刊。

## Q2
* Q:我现在想基于甲状腺癌方面, 做一个研究计划, 最好有一些生物信息分析, 支持研究计划, 并为后续研究提供基础或者线索, 甚至如果能直接得出比较好的结论会更好.

目前，我的想法是：1、收集TCGA甲状腺癌多组学数据库+其他团队发表的甲状腺癌单细胞测序数据；2、将老年甲状腺癌（≥60岁，岁数可以按照数据库样本的分布再议，也可以50岁）挑出，利用聚类或者其他方式将老年甲状腺癌进行分子分型的分类，比如可以分为A亚型、B亚型、C亚型，分析各个亚型间的突变谱、细胞增殖功能、免疫微环境（免疫抑制，免疫激活）、肿瘤细胞粘附能力、淋巴结转移、预后等的差异性；3、根据亚型的特性，比如A亚型偏惰性，那么这类病人就可以暂时保守治疗，其他亚型偏免疫抑制性或者更易转移（内皮细胞粘附能力高等），这类建议积极手术，等等，回答类似的临床问题，为临床手术决策给参考意义（毕竟老年人，开刀风险等等，有些病人可能不开刀定期随访也就是这样）；4、再根据亚型的特点，去决定后续的机制倾向于研究什么，或者对应哪个通路激活就研究哪个机制

请基于以上的内容, 帮我出具详细的可行研究方案和计划

* A:
你的这个思路本质上是在做一件很有潜在临床转化价值的事情：
👉 用 multi-omics + 单细胞数据，在老年甲状腺癌人群中建立“分子分型 → 临床决策”桥梁。

但我要先给你一个关键判断（很重要）：

TCGA-THCA 的“死亡事件数”和“随访时间”是偏弱的，直接做“生存预测+临床决策”可能说服力不足。

👉 所以策略必须调整为：
“机制+分型 + surrogate临床指标（侵袭性/转移/免疫状态）”优先，生存分析作为辅助
这样论文更稳、逻辑更严谨。

⸻

一、整体研究框架（升级版）

🎯 核心科学问题

老年甲状腺癌是否存在分子层面的异质性亚型，这些亚型是否对应不同的：

	•	侵袭性
	•	免疫微环境
	•	转移风险
	•	潜在治疗策略

⸻

🧠 研究逻辑主线（建议这样讲故事）

Step 1

构建 elderly-specific cohort（老年特异队列）

Step 2

做 unsupervised molecular subtyping（无监督分型）

Step 3

定义每个亚型的：
	•	tumor aggressiveness
	•	immune landscape
	•	mutation pattern

Step 4

通过 single-cell数据解释亚型来源机制

Step 5

构建 临床决策导向模型（risk stratification）

⸻

二、数据层设计（非常关键）

1️⃣ Bulk multi-omics（主分析）
	•	TCGA-THCA
	•	RNA-seq
	•	mutation (MAF)
	•	CNV
	•	methylation（可选）
	•	clinical

👉 分组：
	•	≥60岁（建议 primary）
	•	50–60岁（可作为 sensitivity analysis）

⸻

2️⃣ 单细胞数据（机制解释）

来源：Gene Expression Omnibus

建议筛选：
	•	PTC scRNA-seq（tumor + normal + metastasis）
	•	有 immune + stromal annotation 的数据集

👉 用途：
	•	推断亚型来源（tumor cell state vs microenvironment）
	•	验证免疫抑制/激活

⸻

3️⃣ 外部验证（建议加分项）
	•	GEO bulk dataset（验证分型 reproducibility）
	•	或自行收集 cohort（如果你有资源）

⸻

三、核心方法设计（论文价值核心）

⸻

🔬 Aim 1：老年甲状腺癌分子分型

方法：
	•	输入：RNA-seq（TPM/FPKM → log2）
	•	筛选：
	•	高变基因（top 3000–5000）

分型方法（建议组合）：
	•	Consensus Clustering（首选）
	•	NMF（备选）
	•	iCluster（如果做multi-omics integration）

👉 输出：
	•	k=2–4个亚型（不要太多）

⸻

📊 Aim 2：亚型特征解析

1. 突变谱
	•	BRAF / RAS / TERT promoter
	•	TMB
	•	共突变模式（co-mutation）

⸻

2. 肿瘤生物学功能
	•	GSVA / ssGSEA：
	•	proliferation
	•	EMT
	•	angiogenesis
	•	adhesion

👉 特别建议：
	•	adhesion / ECM / integrin（和转移强相关）

⸻

3. 免疫微环境（重点）

方法：
	•	CIBERSORTx
	•	xCell / EPIC

分析：
	•	CD8+ T cells
	•	Treg
	•	macrophage M2
	•	exhaustion markers（PDCD1, LAG3）

👉 定义：
	•	immune-hot
	•	immune-cold
	•	immune-suppressed

⸻

4. 临床关联（关键）

由于生存不足 → 用 surrogate：
	•	tumor stage
	•	lymph node metastasis
	•	extrathyroidal extension
	•	recurrence

👉 用 logistic regression / ML model

⸻

🧬 Aim 3：单细胞机制解析（亮点）

核心问题：

这些分型是由什么驱动？

分析：
	•	细胞类型比例（T cell / CAF / endothelial）
	•	tumor cell state（differentiation vs dedifferentiation）
	•	cell-cell communication（CellChat / NicheNet）

👉 输出：
	•	亚型A：immune active + T cell infiltration
	•	亚型B：CAF enriched + immunosuppressive
	•	亚型C：proliferative tumor dominant

⸻

🧠 Aim 4：临床决策模型（你的核心创新点）

构建：
	•	subtype + clinical features
	•	或 gene signature（LASSO）

输出：
	•	risk score

⸻

最终目标（论文核心图）

👉 一个 decision framework：

亚型	特征	临床建议
A（indolent）	低增殖 + immune active	active surveillance
B（immune suppressed）	Treg/M2高	考虑免疫治疗
C（aggressive）	EMT + adhesion高	建议手术


⸻

四、关键创新点（帮你提升发表层级）

⭐ 创新点1

elderly-specific tumor biology（很少人专门做）

⸻

⭐ 创新点2

multi-omics + single-cell integration

⸻

⭐ 创新点3（最重要）

👉 clinical decision-oriented classification
（不是单纯分型，而是指导是否手术）

⸻

五、潜在风险与规避（非常重要）

❗问题1：TCGA生存数据弱

👉 解决：
	•	用 metastasis / stage 替代

⸻

❗问题2：样本量不足（老年）

👉 解决：
	•	合并50+ vs 60+ 做敏感性分析

⸻

❗问题3：分型不稳定

👉 解决：
	•	bootstrap + external validation

⸻

❗问题4：临床建议“过度推断”

👉 解决（关键）：

不直接说“建议手术”，而是：

	•	“suggest potential stratification for clinical decision-making”

⸻

六、Figure设计（直接给你Nature风格）

Fig1

cohort + clustering

Fig2

mutation landscape

Fig3

immune landscape

Fig4

tumor biology（GSVA）

Fig5

single-cell机制

Fig6

clinical association

Fig7（核心）

👉 decision model / clinical framework

⸻

七、可发表期刊（现实评估）

如果执行到位：

中高水平：
	•	Cancers
	•	Frontiers in Oncology

更高（需要验证数据）：
	•	Clinical Cancer Research
	•	iMeta

⸻

八、我给你的“升级建议”（非常关键）

你当前方案已经很好，但如果你想提高命中率，我建议加一条：

👉 把“免疫逃逸 or ECM/黏附”作为主线机制

原因：
	•	甲状腺癌死亡率低 → 增殖不是关键
	•	转移 + immune escape 才是核心

⸻

九、如果你下一步要继续（我可以帮你做）

我可以帮你直接落地：

1️⃣ 数据下载 + preprocessing pipeline

2️⃣ clustering + 亚型确定代码（R / Python）

3️⃣ immune infiltration + GSVA标准流程

4️⃣ 单细胞整合分析框架

5️⃣ Figure草图（投稿级）

⸻

如果你愿意，我们下一步可以直接进入**“具体分析 pipeline + 代码框架设计”**，甚至可以按你现有服务器环境（你有HPC）直接帮你搭一个可运行的workflow。
