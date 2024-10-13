STEP1安装包

# 设置国内镜像，安装时运行一次即可
options("repos"="https://mirrors.ustc.edu.cn/CRAN/")

if(!"BiocManager"%in%installed.packages()) 
  install.packages("BiocManager",update = F,ask = F)

options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")

# 如果使用R4.0版本需要安装Rtools40
# 下载网站https://cran.r-project.org/bin/windows/Rtools/

# 安装GO.db
if(!"GO.db"%in%installed.packages()) 
  BiocManager::install("GO.db")

# 安装flashClust
if(!"flashClust"%in%installed.packages()) 
  BiocManager::install("flashClust")

# 安装WGCNA
if(!"WGCNA"%in%installed.packages()) 
  BiocManager::install("WGCNA")

# 安装org.Hs.eg.db
if(!"org.Hs.eg.db"%in%installed.packages()) 
  BiocManager::install("org.Hs.eg.db")

# 安装clusterProfiler
if(!"clusterProfiler"%in%installed.packages()) 
  BiocManager::install("clusterProfiler")

# 安装最新版rlang
BiocManager::install("rlang")

# 安装富集分析包
source("https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/GeneAnnotation/installAnRichment.R");
# 上面无法加载成功时，用以下代码本地加载
source("installAnRichment.R")
installAnRichment()

install.packages("anRichmentMethods_0.91-94.tar.gz", repos = NULL, type = "source")
install.packages("anRichment_1.10-1.tar.gz", repos = NULL, type = "source")

# 判断当前目录下是否存在data文件夹
# 存在时，忽略
# 不存在，创建data文件夹储存输入文件和结果文件
if(!dir.exists("data")) dir.create("data")
# 测试figures文件夹是否存在，
# 存在时忽略，不存在时创建
# figures文件夹储存所有的结果图片
if(!dir.exists("figures")) dir.create("figures")


----------------------------------------------------------------------------------------------------------------------------
  
  STEP2数据预处理


rm(list = ls())
library(WGCNA)

# 读取基因表达矩阵数据
fpkm = read.table("data/fpkm.txt", header = T, row.names = 1, check.names = F)
head(fpkm)

### 选取基因方法 ###

## 第一种，通过标准差选择
## 计算每个基因的标准差
fpkm_sd = apply(fpkm,1,sd)#1是对每一行，2是对每一列
## 使用标准差对基因进行降序排序
fpkm_sd_sorted = order(fpkm_sd, decreasing = T)
## 选择前5000个标准差较大的基因
fpkm_num = fpkm_sd_sorted[1:5000]
## 从表达矩阵提取基因
fpkm_filter = fpkm[fpkm_num,]
## 对表达矩阵进行转置
WGCNA_matrix = t(fpkm_filter)#变成行名是样本，列名是基因
## 保存过滤后的数据
save(WGCNA_matrix, file = "data/Step01-fpkm_sd_filter.Rda")

## 第二种，使用绝对中位差选择，推荐使用绝对中位差
WGCNA_matrix = t(fpkm[order(apply(fpkm,1,mad), decreasing = T)[1:5000],])#mad代表绝对中位差
save(WGCNA_matrix, file = "data/Step01-fpkm_mad_filter.Rda")

## 第三种，使用全基因
WGCNA_matrix = t(fpkm)
save(WGCNA_matrix, file = "data/Step01-fpkm_allgene.Rda")


### 去除缺失值较多的基因/样本 ###
rm(list = ls())
# 加载表达矩阵
load(file = "data/Step01-fpkm_mad_filter.Rda")
# 加载WGCNA包
library(WGCNA)
# 判断是否缺失较多的样本和基因
datExpr0 = WGCNA_matrix
gsg = goodSamplesGenes(datExpr0, verbose = 3)

# 是否需要过滤，TRUE表示不需要，FALSE表示需要
gsg$allOK
# 当gsg$allOK为TRUE，以下代码不会运行，为FALSE时，运行以下代码过滤样本
if (!gsg$allOK)
{
  # 打印移除的基因名和样本名
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # 提取保留的基因和样本
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

### 通过样本聚类识别离群样本，去除离群样本 ###
sampleTree = hclust(dist(datExpr0), method = "average");#使用hclust函数进行均值聚类
# 绘制样本聚类图确定离群样本
sizeGrWindow(30,9)
pdf(file = "figures/Step01-sampleClustering.pdf", width = 30, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
# 根据上图判断，需要截取的高度参数h
abline(h = 120, col = "red")#在120的地方画条线
dev.off()

# 去除离群得聚类样本，cutHeight参数要与上述得h参数值一致
clust = cutreeStatic(sampleTree, cutHeight = 120, minSize = 10)
table(clust)
# clust
# 0   1  0就是要去除的，1就是要保存的
# 15 162

# clust 1聚类中包含着我们想要得样本，将其提取出来
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]

# 记录基因和样本数，方便后续可视化
nGenes = ncol(datExpr)#基因数
nSamples = nrow(datExpr)#样本数
save(datExpr, nGenes, nSamples,file = "data/Step01-WGCNA_input.Rda")

----------------------------------------------------------------------------------------------------------------------------
  
  STEP3

# 清空所有变量
rm(list = ls())
# 加载包
library(WGCNA)
# 允许多线程运行
enableWGCNAThreads()
# 加载表达矩阵
load("data/Step01-WGCNA_input.Rda")

### 选择软阈值 ###
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# 进行网络拓扑分析
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)#β=power,就是软阈值

# 可视化结果
sizeGrWindow(9, 5)
pdf(file = "figures/Step02-SoftThreshold.pdf", width = 9, height = 5);

par(mfrow = c(1,2))#一个画板上，画两个图，一行两列
cex1 = 0.9;
# 无尺度网络阈值得选择
plot(sft$fitIndices[,1], 
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",#x轴
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",#y轴
     main = paste("Scale independence"));#标题
text(sft$fitIndices[,1], 
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,
     cex=cex1,
     col="red");

# 用红线标出R^2的参考值
abline(h=0.90,col="red")
# 平均连接度
plot(sft$fitIndices[,1], 
     sft$fitIndices[,5],
     xlab="Soft Threshold (power)",
     ylab="Mean Connectivity", 
     type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], 
     sft$fitIndices[,5], 
     labels=powers, 
     cex=cex1,
     col="red")

dev.off()

# 无尺度网络检验，验证构建的网络是否是无尺度网络
softpower=sft$powerEstimate

ADJ = abs(cor(datExpr,use="p"))^softpower#相关性取绝对值再幂次
k = as.vector(apply(ADJ,2,sum,na.rm=T))#对ADJ的每一列取和，也就是频次

pdf(file = "figures/Step02-scaleFree.pdf",width = 14)
par(mfrow = c(1,2))

hist(k)#直方图
scaleFreePlot(k,main="Check scale free topology")

dev.off()


### 一步构建网络 ###
net = blockwiseModules(datExpr, #处理好的表达矩阵
                       power = sft$powerEstimate,#选择的软阈值
                       TOMType = "unsigned", #拓扑矩阵类型，none表示邻接矩阵聚类，unsigned最常用，构建无方向
                       minModuleSize = 30,#网络模块包含的最少基因数
                       reassignThreshold = 0, #模块间基因重分类的阈值
                       mergeCutHeight = 0.25,#合并相异度低于0.25的模块
                       numericLabels = TRUE, #true，返回模块的数字标签 false返回模块的颜色标签
                       pamRespectsDendro = FALSE,#调用动态剪切树算法识别网络模块后，进行第二次的模块比较，合并相关性高的模块
                       saveTOMs = TRUE,#保存拓扑矩阵
                       saveTOMFileBase = "data/Step02-fpkmTOM", 
                       verbose = 3)#0，不反回任何信息，＞0返回计算过程
# 保存网络构建结果
save(net, file = "data/Step02-One_step_net.Rda")

# 加载网络构建结果
load(file = "data/Step02-One_step_net.Rda")
# 打开绘图窗口
sizeGrWindow(12, 9)
pdf(file = "figures/Step02-moduleCluster.pdf", width = 12, height = 9);
# 将标签转化为颜色
mergedColors = labels2colors(net$colors)
# 绘制聚类和网络模块对应图
plotDendroAndColors(dendro = net$dendrograms[[1]], #hclust函数生成的聚类结果
                    colors = mergedColors[net$blockGenes[[1]]],#基因对应的模块颜色
                    groupLabels = "Module colors",#分组标签
                    dendroLabels = FALSE, #false,不显示聚类图的每个分支名称
                    hang = 0.03,#调整聚类图分支所占的高度
                    addGuide = TRUE, #为聚类图添加辅助线
                    guideHang = 0.05,#辅助线所在高度
                    main = "Gene dendrogram and module colors")
dev.off()

# 加载TOM矩阵
load("data/Step02-fpkmTOM-block.1.RData")
# 网络特征向量
MEs = moduleEigengenes(datExpr, mergedColors)$eigengenes
# 对特征向量排序
MEs = orderMEs(MEs)
# 可视化模块间的相关性
sizeGrWindow(5,7.5);
pdf(file = "figures/Step02-moduleCor.pdf", width = 5, height = 7.5);

par(cex = 0.9)
plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), 
                      marHeatmap = c(3,4,1,2), cex.lab = 0.8, 
                      xLabelsAngle = 90)
dev.off()

## TOMplot

dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = sft$powerEstimate); #1-相关性=相异性

nSelect = 400
# 随机选取400个基因进行可视化，设置seed值，保证结果的可重复性
set.seed(10);#设置随机种子数
select = sample(nGenes, size = nSelect);#从5000个基因选择400个
selectTOM = dissTOM[select, select];#选择这400*400的矩阵
# 对选取的基因进行重新聚类
selectTree = hclust(as.dist(selectTOM), method = "average")#用hclust重新聚类
selectColors = mergedColors[select];#提取相应的颜色模块
# 打开一个绘图窗口
sizeGrWindow(9,9)
pdf(file = "figures/Step02-TOMplot.pdf", width = 9, height = 9);
# 美化图形的设置
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
# 绘制TOM图
TOMplot(plotDiss, #拓扑矩阵，该矩阵记录了每个节点之间的相关性
        selectTree, #基因的聚类结果
        selectColors, #基因对应的模块颜色
        main = "Network heatmap plot, selected genes")
dev.off()
----------------------------------------------------------------------------------------------------------------------------
  
  STEP4

# 清空环境变量
rm(list = ls())
# 加载包
library(WGCNA)
# 加载表达矩阵
load("data/Step01-WGCNA_input.Rda")

# 读入临床信息
clinical = read.table("data/clinical.txt",stringsAsFactors = TRUE, header = T,row.names = 1,na.strings = "",sep = "\t")
# 查看临床信息
head(clinical)
# 对表达矩阵进行预处理
datTraits = as.data.frame(do.call(cbind,lapply(clinical, as.numeric)))
rownames(datTraits) = rownames(clinical)

# 对样本进行聚类
sampleTree2 = hclust(dist(datExpr), method = "average")

# 将临床信息转换为颜色，白色表示低，红色表示高，灰色表示缺失
traitColors = numbers2colors(datTraits, signed = FALSE)

pdf(file = "figures/Step04-Sample_dendrogram_and_trait_heatmap.pdf", width = 24);
# 样本聚类图与样本性状热图
plotDendroAndColors(sampleTree2, 
                    traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")
dev.off()

#### 网络的分析
###### 基因模块与临床信息的关系
# 加载构建的网络
load(file = "data/Step03-Step_by_step_buildnetwork.rda")
# 对模块特征矩阵进行排序
MEs=orderMEs(MEs)
#计算模型特征矩阵和样本信息矩阵的相关度。
moduleTraitCor=cor(MEs, datTraits, use="p")
write.table(file="data/Step04-modPhysiological.cor.xls",moduleTraitCor,sep="\t",quote=F)
moduleTraitPvalue=corPvalueStudent(moduleTraitCor, nSamples)
write.table(file="data/Step04-modPhysiological.p.xls",moduleTraitPvalue,sep="\t",quote=F)

#使用labeledHeatmap()将上述相关矩阵和p值可视化。
pdf(file="figures/Step04-Module_trait_relationships.pdf",width=9,height=7)
textMatrix=paste(signif(moduleTraitCor,2),"\n(",signif(moduleTraitPvalue,1),")",sep="")
dim(textMatrix)=dim(moduleTraitCor)
# 基因模块与临床信息相关性图
labeledHeatmap(Matrix=moduleTraitCor,#模块和表型的相关性矩阵，这个参数最重要，其他可以不变
               xLabels=colnames(datTraits),
               yLabels=names(MEs),
               ySymbols=names(MEs),
               colorLabels=FALSE,
               colors=blueWhiteRed(50),
               textMatrix=textMatrix,
               setStdMargins=FALSE,
               cex.text=0.5,
               cex.lab=0.5,
               zlim=c(-1,1),
               main=paste("Module-trait relationships"))
dev.off()


## 单一模块与某一表型相关性
M_stage = as.data.frame(datTraits$M_stage)
# 分析自己感兴趣的临床信息，此处以M_stage为示例
names(M_stage) = "M_stage"
# 模块对应的颜色
modNames = substring(names(MEs), 3)

# 计算基因模块特征
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

# 对结果进行命名
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

# 计算M分期基因特征显著性
geneTraitSignificance = as.data.frame(cor(datExpr, M_stage, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

# 对结果进行命名
names(geneTraitSignificance) = paste("GS.", names(M_stage), sep="")
names(GSPvalue) = paste("p.GS.", names(M_stage), sep="")

# 设置需要分析的模块名称，此处为brown模块
module = "brown"
# 提取brown模块数据
column = match(module, modNames);
moduleGenes = moduleColors==module;

# 可视化brown模块与M分期的相关性分析结果
sizeGrWindow(7, 7);
pdf(file="figures/Step04-Module_membership_vs_gene_significance.pdf")
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for M Stage",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

## 单一特征与所有模块关联分析
GS = as.numeric(cor(datTraits$M_stage,datExpr, use="p"))
GeneSignificance = abs(GS);
pdf(file="figures/Step04-Gene_significance_for_M_stage_across_module.pdf",width=9,height=5)
plotModuleSignificance(GeneSignificance, moduleColors, ylim=c(0,0.2), main="Gene significance for M stage across module");
dev.off()

## 模块中的hub基因
#### 为每一个模块寻找hub基因
HubGenes <- chooseTopHubInEachModule(datExpr,#WGCNA分析输入的表达矩阵
                                     moduleColors)#模块颜色信息
# 保存hub基因结果
write.table (HubGenes,file = "data/Step04-HubGenes_of_each_module.xls",quote=F,sep='\t',col.names = F)

#### 与某种特征相关的hub基因
NS = networkScreening(datTraits$M_stage,#M分期
                      MEs,#
                      datExpr)#WGCNA分析输入的表达矩阵
# 将结果写入到文件
write.table(NS,file="data/Step04-Genes_for_M_stage.xls",quote=F,sep='\t')

## 模块GO/KEGG分析
# 加载R包
library(anRichment)
library(clusterProfiler)
##### GO分析
# 构建GO背景基因集
GOcollection = buildGOcollection(organism = "human")
geneNames = colnames(datExpr)
# 将基因SYMBOL转换为ENTREZID基因名
geneID = bitr(geneNames,fromType = "SYMBOL", toType = "ENTREZID", 
              OrgDb = "org.Hs.eg.db", drop = FALSE)
# 将基因名对应结果写入文件中
write.table(geneID, file = "data/Step04-geneID_map.xls", sep = "\t", quote = TRUE, row.names = FALSE)
# 进行GO富集分析
GOenr = enrichmentAnalysis(classLabels = moduleColors,#基因所在的模块信息
                           identifiers = geneID$ENTREZID,
                           refCollection = GOcollection,
                           useBackground = "given",
                           threshold = 1e-4,
                           thresholdType = "Bonferroni",
                           getOverlapEntrez = TRUE,
                           getOverlapSymbols = TRUE,
                           ignoreLabels = "grey");
# 提取结果，并写入结果到文件
tab = GOenr$enrichmentTable
names(tab)
write.table(tab, file = "data/Step04-GOEnrichmentTable.xls", sep = "\t", quote = TRUE, row.names = FALSE)

# 提取主要结果，并写入文件
keepCols = c(1, 3, 4, 6, 7, 8, 13)
screenTab = tab[, keepCols]
# 小数位为2位
numCols = c(4, 5, 6)
screenTab[, numCols] = signif(apply(screenTab[, numCols], 2, as.numeric), 2)

# 给结果命名
colnames(screenTab) = c("module", "GOID", "term name", "p-val", "Bonf", "FDR", "size")
rownames(screenTab) = NULL

# 查看结果
head(screenTab)
# 写入文件中
write.table(screenTab, file = "data/Step04-GOEnrichmentTableScreen.xls", sep = "\t", quote = TRUE, row.names = FALSE)

##### KEGG富集分析
# AnRichment没有直接提供KEGG数据的背景集
# 这里使用MSigDBCollection构建C2通路数据集
KEGGCollection = MSigDBCollection("data/msigdb_v7.1.xml", MSDBVersion = "7.1",
                                  organism = "human",
                                  excludeCategories = c("h","C1","C3","C4","C5","C6","C7")) 
# KEGG分析
KEGGenr = enrichmentAnalysis(classLabels = moduleColors,
                             identifiers = geneID$ENTREZID,
                             refCollection = KEGGCollection,
                             useBackground = "given",
                             threshold = 1e-4,
                             thresholdType = "Bonferroni",
                             getOverlapEntrez = TRUE,
                             getOverlapSymbols = TRUE,
                             ignoreLabels = "grey")
# 提取KEGG结果，并写入文件
tab = KEGGenr$enrichmentTable
names(tab)
write.table(tab, file = "data/Step04-KEGGEnrichmentTable.xls", sep = "\t", quote = TRUE, row.names = FALSE)

# 提取主要结果并写入文件
keepCols = c(1, 3, 4, 6, 7, 8, 13)
screenTab = tab[, keepCols]
# 取两位有效数字
numCols = c(4, 5, 6)
screenTab[, numCols] = signif(apply(screenTab[, numCols], 2, as.numeric), 2)

# 对结果表格进行重命名
colnames(screenTab) = c("module", "ID", "term name", "p-val", "Bonf", "FDR", "size")
rownames(screenTab) = NULL
# 查看结果
head(screenTab)
# 写入文件中
write.table(screenTab, file = "data/Step04-KEGGEnrichmentTableScreen.xls", sep = "\t", quote = TRUE, row.names = FALSE)

### 输出cytoscape可视化
# 重新计算TOM，power值设置为前面选择好的
TOM = TOMsimilarityFromExpr(datExpr, power = 7)
# 输出全部网络模块
cyt = exportNetworkToCytoscape(TOM,
                               edgeFile = "data/Step04-CytoscapeInput-edges-all.txt",#基因间的共表达关系
                               nodeFile = "data/Step04-CytoscapeInput-nodes-all.txt",#
                               weighted = TRUE,
                               threshold = 0.1,
                               nodeNames = geneID$SYMBOL,
                               altNodeNames = geneID$ENTREZID,
                               nodeAttr = moduleColors)
# 输出感兴趣网络模块
modules = c("brown", "red")
# 选择上面模块中包含的基因
inModule = is.finite(match(moduleColors, modules))
modGenes = geneID[inModule,]
# 选择指定模块的TOM矩阵
modTOM = TOM[inModule, inModule]

# 输出为Cytoscape软件可识别格式
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("data/Step04-CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("data/Step04-CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.2,
                               nodeNames = modGenes$SYMBOL,
                               altNodeNames = modGenes$ENTREZID,
                               nodeAttr = moduleColors[inModule])

