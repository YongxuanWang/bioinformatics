

####TSMR标准流程####

####初始环境设置####

#文件设置

#GWAS 数据暴露因素和结局简称、标识符、参考基因组标识符、人群类别
EXPOSURE_NAME <- 'BMI'     # 数据名称。建议使用英文简称，如 WBC、BMI 等
EXPOSURE_DATA_ID <- '2'  # 数据标识符。建议使用英文简称，如 UK_Biobank_BMI 等
EXPOSURE_REFER_ID <- '3' # GWAS 数据所用的参考序列版本，如 GRCh37 等
EXPOSURE_POP <- '4'      # 人口学种群标识，如 AFR（非洲）、AMR（美洲）、EAS（东亚）、EUR（欧洲）、SAS（南亚）

OUTCOME_NAME <- '5'
OUTCOME_DATA_ID <- '6'
OUTCOME_REFER_ID <- '7'
OUTCOME_POP <- '8'

# 千人基因组计划中用于 LD 参考的数据集位置
LD_REF <- ''  # 下载 http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz 文件，解压归档。保证该目录下有可用的 .bed、.bim 和 .fam 文件

#创建结果存放根文件目录
root_dir <- "C:/Users/12393/Desktop/MRtest"

# 子文件夹名
topic <- paste0(EXPOSURE_NAME, '(', EXPOSURE_DATA_ID, ')', '→', OUTCOME_NAME, '(', OUTCOME_DATA_ID, ')')
fmt <- '#_FORMATED_DATA'                   # 数据格式化
cor <- '1_correlation_analysis'            # 相关性分析
ld <- '2_linkage_disequilibrium_analysis'  # 连锁不平衡分析
mr_weak <- '3_remove_weak_IV'              # 弱工具变量剔除
mr_con <- '4_remove_confounder'            # 混杂因素剔除
mr <- '5_do_MR'                            # 孟德尔随机化分析
# 子文件夹路径
FMT_DIR <- file.path(root_dir, fmt)
TOPIC_DIR <- file.path(root_dir, topic)
COR_DIR <- file.path(root_dir, topic, cor)
LD_DIR <- file.path(root_dir, topic, ld)
MR_WEAK_DIR <- file.path(root_dir, topic, mr_weak)
MR_CON_DIR <- file.path(root_dir, topic, mr_con)
MR_DIR <- file.path(root_dir, topic, mr)
dirs <- list(FMT_DIR, TOPIC_DIR, COR_DIR, LD_DIR, MR_WEAK_DIR, MR_CON_DIR, MR_DIR)



#文件目录与confounder_SNPs.txt文件初始化
for (dir in dirs) { if (!dir.exists(dir)) { dir.create(dir) } }
confounder_file <- file.path(MR_CON_DIR, '#confounder_SNPs.txt')
if (!file.exists(confounder_file)) { file.create(confounder_file) }

#移除无用的变量，仅保留GLOBAL_VAR
GLOBAL_VAR <- c('EXPOSURE_NAME', 'OUTCOME_NAME',  # 暴露因素和结局的简称
                'EXPOSURE_DATA_ID', 'OUTCOME_DATA_ID',  # GWAS 数据标识符
                'EXPOSURE_REFER_ID', 'OUTCOME_REFER_ID',  # 参考基因组标识符
                'EXPOSURE_POP', 'OUTCOME_POP',  # 人群类别
                'LD_REF',  # 用于 LD 参考的数据集位置
                'FMT_DIR', 'COR_DIR', 'LD_DIR', 'MR_WEAK_DIR', 'MR_CON_DIR', 'MR_DIR',  # 子文件夹路径
                'GLOBAL_VAR')
rm(list = setdiff(ls(), GLOBAL_VAR))
#调包
pkgs <- c(  # 需要安装的包列表
  "VariantAnnotation", "gwasglue", "dplyr", "tidyr",
  "CMplot", "TwoSampleMR", "MendelianRandomization",
  "LDlinkR", "ggplot2", "ggforestplot", "ggfunnel",
  "cowplot", "friendly2MR", "plinkbinr", "FastDownloader",
  "FastTraitR", "MRPRESSO", "SNPlocs.Hsapiens.dbSNP155.GRCh37",
  "SNPlocs.Hsapiens.dbSNP155.GRCh38", "tidyverse",
  'MungeSumstats', 'GenomicFiles'
)
#导入包
lapply(pkgs, function(pkg) {
  suppressMessages(library(pkg, character.only = TRUE, quietly = TRUE))
})
#包初始化结束
rm(list = ls())

####GWAS数据格式化####

#环境设定
setwd(FMT_DIR); cat("当前工作目录：", getwd())

#暴露因素
exposure_file <- "ieu-a-2.vcf"
#结局因素
outcome_file <- "ieu-a-7.vcf" 

#format_gwas_data()函数功能定义
format_gwas_data <- function(file, f_fmt,
                             gws_refer_id, exposure_or_outcome_name, exposure_or_outcome_data_id,
                             chromo = '', position = '', eff_allele = '', non_eff_allele = '', beta_val = '', p_val = '',
                             columns_to_remove = c()) {
  formatted_file_name <- paste0(exposure_or_outcome_name, "_", exposure_or_outcome_data_id, "_formatted.csv.gz")
  if (file.exists(formatted_file_name)) {                                       # 格式化后的文件存在则跳过处理
    return(paste0("文件已存在（", formatted_file_name, "）。跳过格式化……"))
  } else { cat("开始读取数据……\n") }
  if_rename <- FALSE
  if (f_fmt == "tsv") { data <- readr::read_tsv(file); if_rename <- TRUE }       # 预读取 TSV 文件，并计划重命名列
  else if (f_fmt == "csv") { data <- readr::read_csv(file); if_rename <- TRUE }  # 预读取 CSV 文件，并计划重命名列
  else if (f_fmt == "vcf") { data <- file }                                      # VCF 文件不需要预读取
  else { stop(paste0("尚未支持 ", f_fmt, " 格式。数据格式化异常退出！")) }
  if (if_rename) {                                                               # 删除特定列并重命名列
    data <- data %>%
      select(-one_of(columns_to_remove)) %>%  # 删除指定列
      select_if(~!all(is.na(.))) %>%          # 删除所有值均为 NA 的列
      dplyr::rename(
        CHR = chromo,         # SNP 所在染色体号
        BP = position,        # SNP 所在染色体位置
        A1 = eff_allele,      # 效应等位基因
        A2 = non_eff_allele,  # 非效应等位基因
        BETA = beta_val,      # beta 值
        P = p_val             # P 值
      )  # 以上 6 列必须严格声明
    cat("预读取已完成，开始格式化……\n")
  }
  # 删除日志文件（文件名中包含“_log_”）
  for (file in list.files(getwd())) { if (grepl(paste0(exposure_or_outcome_name, "_", exposure_or_outcome_data_id, "_formatted_log_.+\\.txt"), file)) { file.remove(file) } }
  MungeSumstats::format_sumstats(  # 执行格式化
    data,
    ref_genome = gws_refer_id,                            # 用于GWAS的参考基因组标识符
    dbSNP = 155,                                          # dbSNP 版本号（144 或 155）
    ignore_multi_trait = TRUE,                            # 忽略多个 P 值
    strand_ambig_filter = FALSE,                          # 删除具有串不明确基因的 SNP
    bi_allelic_filter = FALSE,                            # 去除非双等位基因的 SNP
    allele_flip_check = FALSE,                            # 对照参考基因组检查并翻转等位基因方向
    indels = FALSE,                                       # 是否包含 indel 变异
    nThread = 16,                                         # 并行的线程数
    save_path = file.path(FMT_DIR, formatted_file_name),  # 输出格式化了的文件
    log_mungesumstats_msgs = TRUE,                        # 保存日志信息
    log_folder = getwd(),                                 # 日志文件夹保存路径
  )
  cat(paste0("格式化已完成：", formatted_file_name, "\n"))
}



#格式化暴露数据
format_gwas_data(file = exposure_file, f_fmt = "vcf",
                 gws_refer_id = EXPOSURE_REFER_ID, exposure_or_outcome_name = EXPOSURE_NAME, exposure_or_outcome_data_id = EXPOSURE_DATA_ID,
                 # chromo = "...", position = "...", eff_allele = "...", non_eff_allele = "...", beta_val = "...", p_val = "..."
)

#格式化结局数据
  format_gwas_data(file = outcome_file, f_fmt = "...",
                   gws_refer_id = OUTCOME_REFER_ID, exposure_or_outcome_name = OUTCOME_NAME, exposure_or_outcome_data_id = OUTCOME_DATA_ID,
                   # chromo = "...", position = "...", eff_allele = "...", non_eff_allele = "...", beta_val = "...", p_val = "...",
  )

rm(list = setdiff(ls(), GLOBAL_VAR))  #移除无用的变量仅保留





####筛选强关联的SNPs####

#环境设定
setwd(COR_DIR); cat("当前工作目录：", getwd())

#暴露因素文件名
exposure_vcf_file <- file.path("ieu-a-2.vcf.gz")

#阈值设定
P_THRESHOLD <- 5e-05  #SNPs筛选的阈值（建议值：5e-8）
output_fmt <- "png"   #支持的格式：jpg、pdf、tiff、png

#读取完成了格式化的数据
exposure_data <- TwoSampleMR::read_exposure_data(  # 注意列名！
  filename = exposure_vcf_file, sep = ',',
  snp_col = 'SNP', chr_col = 'CHR', pos_col = 'BP',
  effect_allele_col = 'A1', other_allele_col = 'A2',
  beta_col = 'BETA', se_col = 'SE', eaf_col = 'FRQ', pval_col = 'P',
  samplesize_col = 'N',
)

#必要时补齐samplesize.exposure列
if (all(is.na(exposure_data$samplesize.exposure))) {
  while (TRUE) {
    num <- readline(prompt = "samplesize.exposure 列缺失，请根据数据来源（网页、MataInfo、文献等）将该列填充为统一数值：")
    if (grepl("^\\d+\\.?\\d*$", num)) {
      exposure_data$samplesize.exposure <- as.numeric(num)
      message("已将 ", num, " 填充至 samplesize.exposure 列")
      break
    } else { stop("输入有误！只支持数值型") }
  }
}

#滤去关联性弱的 SNP 并输出到文件
output_data <- subset(exposure_data, pval.exposure < P_THRESHOLD)
cat(paste("剩余", nrow(output_data), "个 SNP"))
write.csv(output_data, file = 'exposure.pvalue.csv', row.names = FALSE); file.create(paste0("exposure.pvalue.csv - P＜", P_THRESHOLD))

#绘制曼哈顿图：准备数据并绘图输出
exposure_data <- exposure_data[, c("SNP", "chr.exposure", "pos.exposure", "pval.exposure")]
colnames(exposure_data) <- c("SNP", "CHR", "BP", "pvalue")
CMplot(exposure_data,
       plot.type = "m",  # m：线性曼哈顿图
       LOG10 = TRUE, threshold = P_THRESHOLD, threshold.lwd = 3, threshold.lty = 1, signal.cex = 0.2,
       chr.den.col = NULL, cex = 0.2, bin.size = 1e5, ylim = c(0, 50), width = 15, height = 9,
       file.output = TRUE, file = output_fmt, verbose = TRUE
)

rm(list = setdiff(ls(), GLOBAL_VAR))  #移除无用的变量
gc()  #垃圾回收


####去除连锁不平衡效应####

#目录设定
setwd(LD_DIR); cat("当前工作目录：", getwd())

#文件读入
exposure_csv_file <- file.path(COR_DIR, "exposure.pvalue.csv")

#阈值设定
CLUMP_KB <- 500  # 连锁互换的距离（建议值：10000）
CLUMP_R2 <- 0.1  # 连锁互换的 R² 阈值（建议值：0.001）

#读取暴露数据
exposure_data <- read.csv(exposure_csv_file, header = TRUE, sep = ",", check.names = FALSE)

#去除连锁不平衡的 SNP 并保存
unclump_snps <- ieugwasr::ld_clump(dat = dplyr::tibble(rsid = exposure_data$SNP, pval = exposure_data$pval.exposure, id = exposure_data$id.exposure),
                                   clump_kb = CLUMP_KB,
                                   clump_r2 = CLUMP_R2,
                                   plink_bin = get_plink_exe(),
                                   bfile = file.path(LD_REF, EXPOSURE_POP))

exposure_data <- exposure_data %>%
  dplyr::inner_join(unclump_snps, by = c("SNP" = "rsid")) %>%
  dplyr::select(names(.))

#提示句
print(paste("剩余", nrow(exposure_data), "个 SNP"))

#输出数据
write.csv(exposure_data, file = "exposure.LD.csv", row.names = FALSE); file.create(paste0("exposure.LD.csv - ", CLUMP_KB, " kb, ", CLUMP_R2, "R²"))

rm(list = setdiff(ls(), GLOBAL_VAR))  #移除无用的变量



####去除弱工具变量####

#工作目录
setwd(MR_WEAK_DIR); cat("当前工作目录：", getwd())

#写入
exposure_csv_file <- file.path(LD_DIR, "exposure.LD.csv")

#阈值设定
F_THRESHOLD <- 10  # F统计量阈值

#读取连锁不平衡分析结果文件
exposure_data <- read.csv(exposure_csv_file, header = TRUE, sep = ",", check.names = FALSE)

#补齐 R² 和 F 列
exposure_data$R2 <- TwoSampleMR::get_r_from_bsen(exposure_data$beta.exposure, exposure_data$se.exposure, exposure_data$samplesize.exposure)^2
exposure_data$Fval <- (exposure_data$samplesize.exposure - 2) * exposure_data$R2 / (1 - exposure_data$R2)

#过滤保留 F>10 的工具变量
exposure_data <- exposure_data[exposure_data$F > F_THRESHOLD,]
print(paste("剩余", nrow(exposure_data), "个 SNP"))
#写出
write.csv(exposure_data, "exposure.F.csv", row.names = FALSE); file.create(paste0("exposure.F.csv - F＞", F_THRESHOLD))

rm(list = setdiff(ls(), GLOBAL_VAR))  # 移除无用的变量


####去除混杂因素####

#工作目录
setwd(MR_CON_DIR); cat("当前工作目录：", getwd())
#写入
exposure_csv_file <- file.path(MR_WEAK_DIR, "exposure.F.csv")

# ---读取去除了弱工具变量的结果文件
exposure_data <- read.csv(exposure_csv_file, header = TRUE, sep = ",", check.names = FALSE)

# 获取当前的工具变量表型
snp_with_trait <- FastTraitR::look_trait(rsids = exposure_data$SNP, out_file = 'check_SNPs_trait.csv')
snp_with_trait_save <- snp_with_trait %>%
  arrange(trait) %>%
  select(trait) %>%
  distinct()  #工具变量表型去重


#保存到文件
writeLines(snp_with_trait_save$trait, 'check_SNPs_trait.txt') 

#提示语
print(paste("当前筛选到的 SNPs 表型描述，按行分隔地保存到了 check_SNPs_trait.txt 文件～"))

#手动整理混杂因素列表
message(paste0("查看 check_SNPs_trait.txt 文件中的表型是否为 [", EXPOSURE_NAME, " → ", OUTCOME_NAME, "] 的混杂因素，\n将混杂因素保存到 ./4_remove_confounder/#confounder_SNPs.txt 文件！"))
if (file.info("#confounder_SNPs.txt")$size == 0) { stop("请手动整理混杂因素列表文件 #confounder_SNPs.txt！") }


#比较并剔除包含在文本文件中的短语的 SNP
confounders <- readLines("#confounder_SNPs.txt")
snp_with_trait$trait <- tolower(snp_with_trait$trait)  #确保trait列文本均为小写
for (confounder in confounders) {
  snp_with_trait <- snp_with_trait[!grepl(tolower(confounder), snp_with_trait$trait),]
}
snp_with_trait <- dplyr::distinct(snp_with_trait, rsid, .keep_all = FALSE)  #去重
exposure_data <- exposure_data %>%
  dplyr::inner_join(snp_with_trait, by = c("SNP" = "rsid")) %>%
  dplyr::select(names(exposure_data))

#提示语
print(paste("剔除混杂因素后，剩余", nrow(exposure_data), "个 SNP"))
print(paste("剩余", nrow(exposure_data), "个 SNP"))

#保存到文件
write.csv(exposure_data, "exposure.confounder.csv", row.names = FALSE)

rm(list = setdiff(ls(), GLOBAL_VAR))  # 移除无用的变量


####孟德尔随机化分析####

#工作目录
setwd(MR_DIR); cat("当前工作目录：", getwd())

#写入
exposure_csv_file <- file.path(MR_CON_DIR, "exposure.confounder.csv")
outcome_file <- file.path(FMT_DIR, paste0(OUTCOME_NAME, "_", OUTCOME_DATA_ID, "_formatted.csv.gz"))


#读取整理好的暴露数据
exposure_data <- read.csv(exposure_csv_file, header = TRUE, sep = ",", check.names = FALSE)

#读取结局数据
outcome_data <- TwoSampleMR::read_outcome_data(  # 注意列名！
  filename = outcome_file, sep = ',',
  snp_col = 'SNP',
  effect_allele_col = 'A1', other_allele_col = 'A2',
  beta_col = 'BETA', se_col = 'SE', eaf_col = 'FRQ', pval_col = 'P'
)

#筛选结局数据中，与工具变量相交集的部分
outcome_table <- merge(exposure_data, outcome_data, by.x = "SNP", by.y = "SNP")

#输出文件
write.csv(outcome_table[, -(2:ncol(exposure_data))], file = "outcome.csv")

rm(outcome_table)

#将暴露数据和结局数据打标签后整合到一个数据框
exposure_data$exposure <- EXPOSURE_NAME
outcome_data$outcome <- OUTCOME_NAME

#工具变量建立
data <- TwoSampleMR::harmonise_data(exposure_data, outcome_data, action = 1)
#输出工具变量
write.csv(data[data$mr_keep == "TRUE",], file = "SNP.csv", row.names = FALSE)

#MR-PRESSO 方法进行异常值检测，得到偏倚的SNP
mr_presso_result <- run_mr_presso(data)
#输出
write.csv(mr_presso_result[[1]]$
            `MR-PRESSO results`$
            `Outlier Test`, file = "outlier_SNPs.csv")


#----执行MR分析----
mr_result <- mr(data)
#输出
write.csv(generate_odds_ratios(mr_result), file = "MR-Result.csv", row.names = FALSE)

#----异质性检验----
write.csv(mr_heterogeneity(data), file = "heterogeneity.csv", row.names = FALSE)

#----多效性检验----
write.csv(mr_pleiotropy_test(data), file = "pleiotropy.csv", row.names = FALSE)

#可视化
pdf(file = "pic.scatter_plot.pdf", width = 7.5, height = 7); mr_scatter_plot(mr_result, data); dev.off()  # 散点图

res_single <- mr_singlesnp(data)

pdf(file = "pic.forest.pdf", width = 7, height = 6.5); mr_forest_plot(res_single); dev.off()  # 森林图

pdf(file = "pic.funnel_plot.pdf", width = 7, height = 6.5); mr_funnel_plot(singlesnp_results = res_single); dev.off()  # 漏斗图

pdf(file = "pic.leaveoneout.pdf", width = 7, height = 6.5); mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(data)); dev.off()  # 敏感性分析

result_file <- file.path(MR_DIR, "MR-Result.csv")
pFilter <- 1  # 是否根据 P 值过滤？1：不过滤；0.05：过滤
draw_forest_map <- function(inputFile = null, forestFile = null, forestCol = null) {
  # 读取输入文件
  rt <- read.csv(inputFile, header = T, sep = ",", check.names = F)
  row.names(rt) <- rt$method
  rt <- rt[rt$pval < pFilter,]
  method <- rownames(rt)
  or <- sprintf("%.3f", rt$"or")
  orLow <- sprintf("%.3f", rt$"or_lci95")
  orHigh <- sprintf("%.3f", rt$"or_uci95")
  OR <- paste0(or, "(", orLow, "-", orHigh, ")")
  pVal <- ifelse(rt$pval < 0.001, "<0.001", sprintf("%.3f", rt$pval))
  
  #输出图形
  pdf(file = forestFile, width = 7, height = 4.6)
  n <- nrow(rt)
  nRow <- n + 1
  ylim <- c(1, nRow)
  layout(matrix(c(1, 2), nc = 2), width = c(3.5, 2))
  
  # 左侧
  xlim <- c(0, 3)
  par(mar = c(4, 2.5, 2, 1))
  plot(1, xlim = xlim, ylim = ylim, type = "n", axes = F, xlab = "", ylab = "")
  text.cex <- 0.8
  text(0, n:1, method, adj = 0, cex = text.cex)
  text(1.9, n:1, pVal, adj = 1, cex = text.cex); text(1.9, n + 1, 'pvalue', cex = 1, font = 2, adj = 1)
  text(3.1, n:1, OR, adj = 1, cex = text.cex); text(2.7, n + 1, 'OR', cex = 1, font = 2, adj = 1)
  
  # 右侧
  par(mar = c(4, 1, 2, 1), mgp = c(2, 0.5, 0))
  xlim <- c(min(as.numeric(orLow) * 0.975, as.numeric(orHigh) * 0.975, 0.9), max(as.numeric(orLow), as.numeric(orHigh)) * 1.025)
  plot(1, xlim = xlim, ylim = ylim, type = "n", axes = F, ylab = "", xaxs = "i", xlab = "OR")
  arrows(as.numeric(orLow), n:1, as.numeric(orHigh), n:1, angle = 90, code = 3, length = 0.05, col = "darkblue", lwd = 3)
  abline(v = 1, col = "black", lty = 2, lwd = 2)
  boxcolor <- ifelse(as.numeric(or) > 1, forestCol, forestCol)
  points(as.numeric(or), n:1, pch = 15, col = boxcolor, cex = 2)
  axis(1)
  dev.off()
}
draw_forest_map(inputFile = result_file, forestFile = "forest_map.pdf", forestCol = "red")

rm(list = setdiff(ls(), GLOBAL_VAR))  # 移除无用的变量

#移除无用的变量
rm(list = setdiff(ls(), GLOBAL_VAR))  
gc()  # 垃圾回收






####双样本孟德尔随机化简略版####
#先找到暴露因素和结局的ID号
#ckd：finn-b-N14_CHRONKIDNEYDIS   nicm：finn-b-I9_NONISCHCARDMYOP
#安装TSMR包
install.packages("devtools")
devtools::install_github("MRCIEU/TwoSampleMR",force = TRUE)

#运行TSMR
library(TwoSampleMR)
#提取暴露因素的工作变量
exposure_input <- extract_instruments(outcomes = 'finn-b-N14_CHRONKIDNEYDIS')

#解释extract_instruments函数
extract_instruments(outcomes,p1 = 5e-08,clump = TRUE,
                    p2 = 5e-08,r2 = 0.001,kb = 10000,access_token = ieugwasr::check_access_token(),
                    force_server = FALSE)

#如果想要调整clump参数
exposure_input <- extract_instruments(outcomes = 'ebi-a-GCST90018866',
                                      clump = TRUE, r2 = 0.01,
                                      kb = 5000, access_token = NULL)

#如果想要调整Pֵ值
exposure_input <- extract_instruments(outcomes = 'ebi-a-GCST008311',
                                      p1 = 5e-06,
                                      clump = TRUE, r2 = 0.001,
                                      kb = 10000, access_token = NULL)


#提取工作变量在结局中的信息
outcome_dat_input <- extract_outcome_data(snps = exposure_input$SNP, outcomes = 'ebi-a-GCST90018902')

#将暴露和结局数据进行合并，产生用于进行MR分析的数据
#第一种代码
dat <- harmonise_data(exposure_input, outcome_dat_input)
#第二种代码
dat <- harmonise_data(
  exposure_dat=bmi_exp_dat,
  outcome_dat=chd_out_dat,
  action= 2
)

#mr分析的主要结果：默认用5种方法进行MR分析
res <- mr(dat)
res

#b值转换成OR值
MRresult <-generate_odds_ratios(res)
MRresult

#如果MR分析中限定方法，如只用mr_egger和mr_ivw
mr(dat, method_list = c("mr_egger_regression", "mr_ivw"))

#ʹ使用随机效应模型
RE <-mr(dat,method_list=c('mr_ivw_mre'))
REOR <-generate_odds_ratios(RE)

#固定效应模型
FE <-mr(dat,method_list=c('mr_ivw_fe'))
FEOR <-generate_odds_ratios(FE)

#离群值检验
#安装MRPRESSO包
devtools::install_github("rondolab/MR-PRESSO",force = TRUE)

#安装完运行
library(MRPRESSO)
mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
          OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = dat, NbDistribution = 1000,  
          SignifThreshold = 0.05)

#单个snp分析（三种方法）
#默认wald比值
res_single <- mr_singlesnp(dat)
ORR <-generate_odds_ratios(res_single)

#数据敏感性分析
#异质性检验
heterogeneity <- mr_heterogeneity(dat)
heterogeneity
#多效性检验
Horizontal_pleiotropy <- mr_pleiotropy_test(dat)
Horizontal_pleiotropy
#留一法检验
single <- mr_leaveoneout(dat)
mr_leaveoneout_plot(single)
#散点图
mr_scatter_plot(res,dat)
#森林图
res_single <- mr_singlesnp(dat)
mr_forest_plot(res_single)
#漏斗图
mr_funnel_plot(res_single)


#保存数据
#安装xlsx
install.packages('xlsx')
#运行xlsx
library(xlsx)
write.xlsx(MRresult, "D:MRresult.xlsx")

#本地为暴露文件，结局为CHD的ID：
#运行TSMR
library(TwoSampleMR)
#读取暴露本地数据
exp_dat <- read_exposure_data(
  filename = 'Blood selenium.csv',
  sep= ",",
  snp_col = "SNP",
  beta_col = "Beta",
  se_col = "SE",
  effect_allele_col ="EA",
  other_allele_col = "NEA",
  eaf_col = "EAF",
  pval_col = "P"
)
exp_dat$exposure <- "Blood selenium"
#读取工具变量在结局当中的信息
CHD_out <- extract_outcome_data(
  snps=exp_dat$SNP,
  outcomes='ieu-a-7',
  proxies = FALSE,
  maf_threshold = 0.01,
  access_token = NULL
)
mydata <- harmonise_data(
  exposure_dat=exp_dat,
  outcome_dat=CHD_out,
  action= 2
)
res <- mr(mydata)
res
OR <-generate_odds_ratios(res)
OR
#异质性分析
het <- mr_heterogeneity(mydata)
het
#多效性检验
pleio <- mr_pleiotropy_test(mydata)
pleio
#留一法
single <- mr_leaveoneout(mydata)
mr_leaveoneout_plot(single)
#散点图
mr_scatter_plot(res,mydata)
#森林图
res_single <- mr_singlesnp(mydata)
mr_forest_plot(res_single)
#漏斗图
mr_funnel_plot(res_single)

#结局为本地文件，怎么读取
#女性肾细胞癌GWAS数据读取
library(data.table)
t2d <-fread('RCC_Females.txt',header=T)
#看数据前六列
head(t2d)

#当将结局变量转换成我们可以使用的TSMR中结局得到形式，应用函数format_data

pFilter=1       #pvalue过滤条件
setwd("G:\\MR\\LAM-pneumothorax")     #设置工作目录

#定义森林图函数
bioForest=function(inputFile=null, forestFile=null, forestCol=null){
  #读取输入文件
  rt=read.csv(inputFile, header=T, sep=",", check.names=F)
  row.names(rt)=rt$method
  rt=rt[rt$pval<pFilter,]
  method <- rownames(rt)
  or <- sprintf("%.3f",rt$"or")
  orLow  <- sprintf("%.3f",rt$"or_lci95")
  orHigh <- sprintf("%.3f",rt$"or_uci95")
  OR <- paste0(or,"(",orLow,"-",orHigh,")")
  pVal <- ifelse(rt$pval<0.001, "<0.001", sprintf("%.3f", rt$pval))
  
  #输出图形
  pdf(file=forestFile, width=7, height=4.6)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3.5,2))
  
  #绘制森林图左边的信息
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,method,adj=0,cex=text.cex)
  text(1.9,n:1,pVal,adj=1,cex=text.cex);text(1.9,n+1,'pvalue',cex=1,font=2,adj=1)
  text(3.1,n:1,OR,adj=1,cex=text.cex);text(2.7,n+1,'OR',cex=1,font=2,adj=1)
  
  #绘制右边的森林图
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(min(as.numeric(orLow)*0.975,as.numeric(orHigh)*0.975,0.9),max(as.numeric(orLow),as.numeric(orHigh))*1.025)
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="OR")
  arrows(as.numeric(orLow),n:1,as.numeric(orHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=3)
  abline(v=1, col="black", lty=2, lwd=2)
  boxcolor = ifelse(as.numeric(or)>1, forestCol, forestCol)
  points(as.numeric(or), n:1, pch = 15, col = boxcolor, cex=2)
  axis(1)
  dev.off()
}

#调用函数，绘制森林图
bioForest(inputFile="MRresult.csv", forestFile="forest.pdf", forestCol="red")





#### 中介MR分析 ####

library(TwoSampleMR)
library(ggplot2)
library(foreach)
source("BWMR_updated.R")


setwd("C:/Users/12393/Desktop/新建文件夹")



#使用循环遍历来完成MR分析
#foreach()可以并行计算
#本段代码使用一个列表存储大量暴露数据id，将id对应暴露数据整合入clump文件内
#逐个读取暴露并对同一个结局数据进行MR分析



##1、批量筛选与结局相关的A因子

#1 
iddf=read.table("lipidomid.txt",header =T,sep = "\t")#读取粗列表
dxw=as.vector(iddf$id)
result <- data.frame()

foreach(i=dxw, .errorhandling = "pass") %do%{
  
  #读取曝露数据
  expo_rt <- read_exposure_data(
    filename = paste0("clump/", i, ".txt"),
    sep = "\t",
    snp_col = "rsid",
    beta_col = "beta",
    se_col = "standard_error",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "effect_allele_frequency",
    pval_col = "pval",
    samplesize_col = "n"
  )
  
  #读取结局数据
  outc_rt <- read_outcome_data(
    snps = expo_rt$SNP,
    filename = "myoutcome.gz",
    sep = "\t",
    snp_col = "rsids",
    beta_col = "beta",
    se_col = "sebeta",
    effect_allele_col = "alt",
    other_allele_col = "ref",
    eaf_col = "af_alt",
    pval_col = "pval"
  )
  
  # 协调处理曝露数据和结果数据
  harm_rt <- harmonise_data(
    exposure_dat = expo_rt,
    outcome_dat = outc_rt,
    action = 2
  )
  
  # 计算R平方（R2）
  harm_rt$R2 <- (2 * (harm_rt$beta.exposure^2) * harm_rt$eaf.exposure * (1 - harm_rt$eaf.exposure)) /
    (2 * (harm_rt$beta.exposure^2) * harm_rt$eaf.exposure * (1 - harm_rt$eaf.exposure) +
       2 * harm_rt$samplesize.exposure * harm_rt$eaf.exposure * (1 - harm_rt$eaf.exposure) * harm_rt$se.exposure^2)
  
  # 计算F统计量
  harm_rt$f <- harm_rt$R2 * (harm_rt$samplesize.exposure - 2) / (1 - harm_rt$R2)
  
  # 计算F统计量的均值
  harm_rt$meanf <- mean(harm_rt$f)
  
  # 筛选出F统计量大于10的数据
  harm_rt <- harm_rt[harm_rt$f > 10, ]
  
  # 执行MR分析
  mr_result <- mr(harm_rt)
  
  # 生成OR结果
  result_or <- generate_odds_ratios(mr_result)
  
  
  # 创建以i命名的文件夹
  filename <- paste0("final_result/", i)
  if (!file.exists(filename)) {
    dir.create(filename)
  }
  
  # 写入文件
  write.table(harm_rt, file = paste0(filename, "/harmonise.txt"), row.names = FALSE, sep = "\t", quote = FALSE)
  write.table(result_or[, 5:ncol(result_or)], file = paste0(filename, "/OR.txt"), row.names = FALSE, sep = "\t", quote = FALSE)
  
  # 如果result_or中第三个元素p-value小于0.05，则将相关信息添加到result中
  if (result_or$pval[3] < 0.05) {
    result <- rbind(result, cbind(id = i, pvalue = result_or$pval[3]))
  }
}

write.table(result,"result.txt",sep = "\t",quote = F,row.names = F)

#2
iddf=read.table("result.txt",header =T,sep = "\t")#读取筛选后的列表
dxw=as.vector(iddf$id)
result1=data.frame()


foreach(i=dxw, .errorhandling = "pass") %do%{
  
  #读取曝露数据
  expo_rt <- read_exposure_data(
    filename = paste0("clump/", i, ".txt"),
    sep = "\t",
    snp_col = "rsid",
    beta_col = "beta",
    se_col = "standard_error",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "effect_allele_frequency",
    pval_col = "pval",
    samplesize_col = "n"
  )
  
  #读取结局数据
  outc_rt <- read_outcome_data(
    snps = expo_rt$SNP,
    filename = "myoutcome.gz",
    sep = "\t",
    snp_col = "rsids",
    beta_col = "beta",
    se_col = "sebeta",
    effect_allele_col = "alt",
    other_allele_col = "ref",
    eaf_col = "af_alt",
    pval_col = "pval"
  )
  
  # 协调处理曝露数据和结果数据
  harm_rt <- harmonise_data(
    exposure_dat = expo_rt,
    outcome_dat = outc_rt,
    action = 2
  )
  
  # 计算R平方（R2）
  harm_rt$R2 <- (2 * (harm_rt$beta.exposure^2) * harm_rt$eaf.exposure * (1 - harm_rt$eaf.exposure)) /
    (2 * (harm_rt$beta.exposure^2) * harm_rt$eaf.exposure * (1 - harm_rt$eaf.exposure) +
       2 * harm_rt$samplesize.exposure * harm_rt$eaf.exposure * (1 - harm_rt$eaf.exposure) * harm_rt$se.exposure^2)
  
  # 计算F统计量
  harm_rt$f <- harm_rt$R2 * (harm_rt$samplesize.exposure - 2) / (1 - harm_rt$R2)
  
  # 计算F统计量的均值
  harm_rt$meanf <- mean(harm_rt$f)
  
  # 筛选出F统计量大于10的数据
  harm_rt <- harm_rt[harm_rt$f > 10, ]
  
  # 执行MR分析
  mr_result <- mr(harm_rt)
  
  # 生成OR结果
  result_or <- generate_odds_ratios(mr_result)
  
  
  # 创建以i命名的文件夹
  filename <- paste0("final_result/", i)
  if (!file.exists(filename)) {
    dir.create(filename)
  }
  
  # 写入文件
  write.table(harm_rt, file = paste0(filename, "/harmonise1.txt"), row.names = FALSE, sep = "\t", quote = FALSE)
  write.table(result_or[, 5:ncol(result_or)], file = paste0(filename, "/OR.txt"), row.names = FALSE, sep = "\t", quote = FALSE)
  harmrt1<-harm_rt
  harmrt2<-harm_rt
  
  
  # 多效性检验
  pleiotropy <- mr_pleiotropy_test(harmrt1)
  write.table(pleiotropy, file = paste0(filename, "/pleiotropy.txt"), sep = "\t", quote = FALSE)
  
  # 异质性检验
  heterogeneity <- mr_heterogeneity(harmrt1)
  write.table(heterogeneity, file = paste0(filename, "/heterogeneity.txt"), sep = "\t", quote = FALSE)
  
  # 散点图绘制
  p1 <- mr_scatter_plot(mr_result, harmrt1)
  ggsave(p1[[1]], file = paste0(filename, "/scatter.pdf"), width = 8, height = 8)
  
  # 单个SNP分析
  singlesnp_res <- mr_singlesnp(harmrt1)
  singlesnpOR <- generate_odds_ratios(singlesnp_res)
  write.table(singlesnpOR, file = paste0(filename, "/singlesnpOR.txt"), row.names = FALSE, sep = "\t", quote = FALSE)
  p4 <- mr_funnel_plot(singlesnp_res)
  ggsave(p2[[1]], file = paste0(filename, "/funnelplot.pdf"), width = 8, height = 8)
  
  # 森林图绘制
  p2 <- mr_forest_plot(singlesnp_res)
  ggsave(p3[[1]], file = paste0(filename, "/forest.pdf"), width = 8, height = 8)
  
  # 逐个剔除检验
  sen_res <- mr_leaveoneout(harmrt1)
  p3 <- mr_leaveoneout_plot(sen_res)
  ggsave(p4[[1]], file = paste0(filename, "/sensitivity-analysis.pdf"), width = 8, height = 8)
  
  
  # MR PRESSO 分析
  presso <- run_mr_presso(harmrt1, NbDistribution = 1000)
  capture.output(presso, file = paste0(filename, "/presso.txt"))
  
  
  # 贝叶斯加权分析
  myBWMR <- BWMR(gammahat = harmrt2$beta.exposure,
                 Gammahat = harmrt2$beta.outcome,
                 sigmaX = harmrt2$se.exposure,
                 sigmaY = harmrt2$se.outcome) 
  beta=myBWMR[["beta"]]
  lci95 <- myBWMR[["beta"]]-1.96*myBWMR[["se_beta"]]
  uci95 <- myBWMR[["beta"]]+1.96*myBWMR[["se_beta"]]
  or <- exp(myBWMR[["beta"]])
  or_lci95 <- exp(lci95)
  or_uci95 <- exp(uci95)
  pval<-myBWMR[["P_value"]]
  result1=rbind(result1,cbind(id=i,method="BWMR",beta,lci95,uci95,or,or_lci95,or_uci95,pval))
  
}

write.table(result1,"BWMR_OR.txt",row.names = F,sep = "\t",quote = F)


##2、批量筛选与A因子相关的B因子
#计算A因子对B因子的重要度
#1
i="GCST90277263"#选取自己的兴趣A因子

dxw=read.table("id91.txt",header = T,sep = "\t")#读取B因子列表
dxwid=as.vector(dxw$id)

#dxwid=dxwid[1:10]

result=data.frame()
foreach(j=dxwid, .errorhandling = "pass") %do%{
  expo_rt=read_exposure_data(file = paste0("clump/",i,".txt"),
                             sep = "\t",
                             snp_col = "rsid",
                             beta_col = "beta",
                             se_col = "standard_error",
                             effect_allele_col = "effect_allele",
                             other_allele_col = "other_allele",
                             eaf_col = "effect_allele_frequency",
                             pval_col = "pval",
                             samplesize_col = "n")
  

  outc_rt <- read_outcome_data(
    snps = expo_rt$SNP,
    filename = paste0(j,".tsv.gz"),
    sep = "\t",
    snp_col = "rsid",
    beta_col = "beta",
    se_col = "standard_error",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "effect_allele_frequency",
    pval_col = "p_value",
    samplesize_col = "n"
  )
  
  
  harm_rt <- harmonise_data(exposure_dat =  expo_rt,outcome_dat = outc_rt,action=2)
  harm_rt$R2 <- (2 * (harm_rt$beta.exposure^2) * harm_rt$eaf.exposure * (1 - harm_rt$eaf.exposure)) /
    (2 * (harm_rt$beta.exposure^2) * harm_rt$eaf.exposure * (1 - harm_rt$eaf.exposure) +
       2 * harm_rt$samplesize.exposure*harm_rt$eaf.exposure * (1 - harm_rt$eaf.exposure) * harm_rt$se.exposure^2)
  harm_rt$f <- harm_rt$R2 * (harm_rt$samplesize.exposure - 2) / (1 - harm_rt$R2)
  harm_rt$meanf<- mean( harm_rt$f)
  harm_rt<-harm_rt[harm_rt$f>10,]
  mr_result<- mr(harm_rt)
  result_or=generate_odds_ratios(mr_result)
  if(result_or$pval[3]<0.05){
    result=rbind(result,cbind(lip=i,inf=j,pvalue=result_or$pval[3]))
    
  }
}
write.table(result,"firstselect.txt",sep = "\t",quote = F,row.names = F)


#2
iddf=read.table("firstselect.txt",header = T,sep = "\t")
iddf=as.vector(iddf$inf)

foreach(j=iddf, .errorhandling = "pass") %do%{
  expo_rt=read_exposure_data(file = paste0("clump/",i,".txt"),
                             sep = "\t",
                             snp_col = "rsid",
                             beta_col = "beta",
                             se_col = "standard_error",
                             effect_allele_col = "effect_allele",
                             other_allele_col = "other_allele",
                             eaf_col = "effect_allele_frequency",
                             pval_col = "pval",samplesize_col = "n")

  outc_rt <- read_outcome_data(
    snps = expo_rt$SNP,
    filename = paste0(j,".tsv.gz"),
    sep = "\t",
    snp_col = "rsid",
    beta_col = "beta",
    se_col = "standard_error",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "effect_allele_frequency",
    pval_col = "p_value",
    samplesize_col = "n"
  )
  harm_rt <- harmonise_data(
    exposure_dat =  expo_rt, 
    outcome_dat = outc_rt,action=2)
  harm_rt$R2 <- (2 * (harm_rt$beta.exposure^2) * harm_rt$eaf.exposure * (1 - harm_rt$eaf.exposure)) /
    (2 * (harm_rt$beta.exposure^2) * harm_rt$eaf.exposure * (1 - harm_rt$eaf.exposure) +
       2 * harm_rt$samplesize.exposure*harm_rt$eaf.exposure * (1 - harm_rt$eaf.exposure) * harm_rt$se.exposure^2)
  harm_rt$f <- harm_rt$R2 * (harm_rt$samplesize.exposure - 2) / (1 - harm_rt$R2)
  harm_rt$meanf<- mean( harm_rt$f)
  harm_rt<-harm_rt[harm_rt$f>10,]
  mr_result<- mr(harm_rt)
  result_or=generate_odds_ratios(mr_result)
  newid=paste0(i,"_",j)
  filename=paste0("final_result/",newid)
  dir.create(filename)
  write.table(harm_rt, file =paste0(filename,"/harmonise2.txt"),row.names = F,sep = "\t",quote = F)
  write.table(result_or[,5:ncol(result_or)],file =paste0(filename,"/OR.txt"),row.names = F,sep = "\t",quote = F)
  pleiotropy=mr_pleiotropy_test(harm_rt)
  write.table(pleiotropy,file = paste0(filename,"/pleiotropy.txt"),sep = "\t",quote = F)
  heterogeneity=mr_heterogeneity(harm_rt)
  write.table(heterogeneity,file = paste0(filename,"/heterogeneity.txt"),sep = "\t",quote = F)
  p1 <- mr_scatter_plot(mr_result, harm_rt)
  ggsave(p1[[1]], file=paste0(filename,"/scatter.pdf"), width=8, height=8)
  singlesnp_res<- mr_singlesnp(harm_rt)
  singlesnpOR=generate_odds_ratios(singlesnp_res)
  write.table(singlesnpOR,file=paste0(filename,"/singlesnpOR.txt"),row.names = F,sep = "\t",quote = F)
  p2 <- mr_forest_plot(singlesnp_res)
  ggsave(p2[[1]], file=paste0(filename,"/forest.pdf"), width=8, height=8)
  sen_res<- mr_leaveoneout(harm_rt)
  p3 <- mr_leaveoneout_plot(sen_res)
  ggsave(p3[[1]], file=paste0(filename,"/sensitivity-analysis.pdf"), width=8, height=8)
  res_single <- mr_singlesnp(harm_rt)
  p4 <- mr_funnel_plot(singlesnp_res)
  ggsave(p4[[1]], file=paste0(filename,"/funnelplot.pdf"), width=8, height=8)
  presso=run_mr_presso(harm_rt,NbDistribution = 1000)
  capture.output(presso,file = paste0(filename,"/presso.txt"))
}



##3、批量筛选与结局相关的B因子

iddf=read.table("firstselect.txt",header =T,sep = "\t")

dxw=as.vector(iddf$inf)
result=data.frame()

foreach(i=dxw, .errorhandling = "pass") %do%{
  expo_rt<- read_exposure_data(
    filename = paste0("clump5/",i,".txt"),
    sep = "\t",
    snp_col = "rsid",
    beta_col = "beta",
    se_col = "standard_error",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "effect_allele_frequency",
    pval_col = "p_value",
    samplesize_col = "n"
  )

  outc_rt <- read_outcome_data(
    snps = expo_rt$SNP,
    filename = "myoutcome.gz",
    sep = "\t",
    snp_col = "rsids",
    beta_col = "beta",
    se_col = "sebeta",
    effect_allele_col = "alt",
    other_allele_col = "ref",
    eaf_col = "af_alt",
    pval_col = "pval")
  
  
  
  harm_rt <- harmonise_data(
    exposure_dat =  expo_rt, 
    outcome_dat = outc_rt,action=2)
  harm_rt$R2 <- (2 * (harm_rt$beta.exposure^2) * harm_rt$eaf.exposure * (1 - harm_rt$eaf.exposure)) /
    (2 * (harm_rt$beta.exposure^2) * harm_rt$eaf.exposure * (1 - harm_rt$eaf.exposure) +
       2 * harm_rt$samplesize.exposure*harm_rt$eaf.exposure * (1 - harm_rt$eaf.exposure) * harm_rt$se.exposure^2)
  harm_rt$f <- harm_rt$R2 * (harm_rt$samplesize.exposure - 2) / (1 - harm_rt$R2)
  harm_rt$meanf<- mean( harm_rt$f)
  harm_rt<-harm_rt[harm_rt$f>10,]
  
  mr_result<- mr(harm_rt)
  result_or=generate_odds_ratios(mr_result) 
  if(result_or$pval[3]<0.05){
    dir.create(i) 
    filename=i
    write.table(harm_rt, file =paste0(filename,"/harmonise3.txt"),row.names = F,sep = "\t",quote = F)
    write.table(result_or[,5:ncol(result_or)],file =paste0(filename,"/OR.txt"),row.names = F,sep = "\t",quote = F)
    pleiotropy=mr_pleiotropy_test(harm_rt)
    write.table(pleiotropy,file = paste0(filename,"/pleiotropy.txt"),sep = "\t",quote = F)
    heterogeneity=mr_heterogeneity(harm_rt)
    write.table(heterogeneity,file = paste0(filename,"/heterogeneity.txt"),sep = "\t",quote = F)

    p1 <- mr_scatter_plot(mr_result, harm_rt)
    ggsave(p1[[1]], file=paste0(filename,"/scatter.pdf"), width=8, height=8)

    singlesnp_res<- mr_singlesnp(harm_rt)
    singlesnpOR=generate_odds_ratios(singlesnp_res)
    write.table(singlesnpOR,file=paste0(filename,"/singlesnpOR.txt"),row.names = F,sep = "\t",quote = F)
    p2 <- mr_forest_plot(singlesnp_res)
    ggsave(p2[[1]], file=paste0(filename,"/forest.pdf"), width=8, height=8)
    sen_res<- mr_leaveoneout(harm_rt)
    p3 <- mr_leaveoneout_plot(sen_res)
    ggsave(p3[[1]], file=paste0(filename,"/sensitivity-analysis.pdf"), width=8, height=8)
    res_single <- mr_singlesnp(harm_rt)
    p4 <- mr_funnel_plot(singlesnp_res)
    ggsave(p4[[1]], file=paste0(filename,"/funnelplot.pdf"), width=8, height=8)
    presso=run_mr_presso(harm_rt,NbDistribution = 1000)
    capture.output(presso,file = paste0(filename,"/presso.txt"))
  }
}


##4、中介效应计算


harm_rt=read.table("harmonise1.txt",sep = "\t",header = T)
mr_result<- mr(harm_rt)
View(mr_result)
result_or=generate_odds_ratios(mr_result)
write.table(result_or[,5:ncol(result_or)],"OR1.txt",row.names = F,sep = "\t",quote = F)

# step 2 outcome to lipidom（反向验证）
#直接看前面的分析


#step 3 lipidom to inf
harm_rt3=read.table("harmonise2.txt",sep = "\t",header = T)
mr_result3<- mr(harm_rt3)
View(mr_result3)
result_or3=generate_odds_ratios(mr_result3)
write.table(result_or3[,5:ncol(result_or3)],"OR3.txt",row.names = F,sep = "\t",quote = F)


#step 4 inf to outcome
harm_rt4=read.table("harmonise3.txt",sep = "\t",header = T)
mr_result4<- mr(harm_rt4)
View(mr_result4)
result_or4=generate_odds_ratios(mr_result4)
write.table(result_or4[,5:ncol(result_or4)],"OR4.txt",row.names = F,sep = "\t",quote = F)


#####mediation effect
beta_all=mr_result[3,"b"]

beta1=mr_result3[3,"b"]

beta2=mr_result4[3,"b"]

#中介效应（间接效应）
beta12=beta1*beta2

#直接效应
beta_dir=beta_all-beta12

se=sqrt(mr_result3[3,"b"]^2*mr_result3[3,"se"]^2+mr_result4[3,"b"]^2*mr_result4[3,"se"]^2)

#Z
Z=beta12/se

#p
P=2*pnorm(q=abs(Z), lower.tail=FALSE)

#####95%可信区间
lci=beta12-1.96*se

uci=beta12+1.96*se

#中介效应所占的比例
beta12_p=beta12/beta_all
lci_p=lci/beta_all
uci_p=uci/beta_all


####eQTL数据提取####
library(data.table)
all_dfd=fread("Whole_Blood.v8.EUR.signif_pairs.txt.gz")#加载数据
head(all_dfd)
#单基因数据提取
mygene="ATAD3B"
myeqtl=all_dfd[all_dfd$external_gene_name==mygene,]
fwrite(myeqtl,file = paste0(mygene,".txt"),sep = "\t",quote = F)
#多基因数据提取
genelist=read.table("mygenelist.txt",header = F,sep = "\t")
genelist1=as.vector(genelist[,1])
foreach(i=genelist1, .errorhandling = "pass") %do%{
  
  myeqtl=all_dfd[all_dfd$external_gene_name==i,]
  
  fwrite(myeqtl,file = paste0(i,".txt"),sep = "\t",quote = F)
}

#r包本地运行加载
library(data.table)
library(foreach)
source("ld_clump.R")
source("ld_matrix.R")
source("afl2.r")
source("api.R")
source("backwards.R")
source("query.R")
source("utils-pipe.R")
source("variants.R")
source("zzz.R")

#基因数据预处理
genelist=read.table("mygenelist.txt",header = F,sep = "\t")
genelist1=as.vector(genelist[,1])

foreach(i=genelist1, .errorhandling = "pass") %do%{
  
  expo_rt=fread(file = paste0(i,".txt"), header = T)
  
  expo_rt=expo_rt[expo_rt$p<1e-5,]

  expo_rt2=expo_rt[,c("rsids2","p")]
  colnames(expo_rt2)=c("rsid", "pval")
  
  #去除LD
  clumdf <- ld_clump_local(dat = expo_rt2, clump_kb = 100, clump_r2 = 0.1,clump_p=1,
                           bfile ="C:\\Users\\77632\\Desktop\\geneMR\\3_step3\\1\\data_maf0\\data_maf0.01_rs_ref", 
                           plink_bin = "C:\\Users\\77632\\Desktop\\geneMR\\3_step3\\1\\plink_win64_20231018\\plink.exe")
  
  #检查
  expo_rt3=expo_rt[which(expo_rt$rsids2%in%clumdf$rsid),]
  
  write.table(expo_rt3,file = paste0(i,"_expo_rt.txt"),row.names = F,sep = "\t",quote = F)

}

#读取暴露数据
expo_rt<- read_exposure_data(
  filename = "_expo_rt.txt",
  sep = "\t",
  snp_col = "rsids2",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "eaf",
  pval_col = "p",
  samplesize_col = "n")
#读取结局数据
outc_rt <- read_outcome_data(
  snps = expo_rt$SNP,
  filename = "myoutcome.gz",
  sep = "\t",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt",
  pval_col = "pval")

#单基因流程
harm_rt <- harmonise_data(
  exposure_dat =  expo_rt, 
  outcome_dat = outc_rt,action=2)

harm_rt$R2 <- (2 * (harm_rt$beta.exposure^2) * harm_rt$eaf.exposure * (1 - harm_rt$eaf.exposure)) /
  (2 * (harm_rt$beta.exposure^2) * harm_rt$eaf.exposure * (1 - harm_rt$eaf.exposure) +
     2 * harm_rt$samplesize.exposure*harm_rt$eaf.exposure * (1 - harm_rt$eaf.exposure) * harm_rt$se.exposure^2)

harm_rt$f <- harm_rt$R2 * (harm_rt$samplesize.exposure - 2) / (1 - harm_rt$R2)

harm_rt$meanf<- mean( harm_rt$f)

harm_rt<-harm_rt[harm_rt$f>10,]

#mr
mr_result<- mr(harm_rt)
result_or=generate_odds_ratios(mr_result) 
write.table(harm_rt, file ="harmonise.txt",row.names = F,sep = "\t",quote = F)
write.table(result_or[,5:ncol(result_or)],file ="OR.txt",row.names = F,sep = "\t",quote = F)
pleiotropy=mr_pleiotropy_test(harm_rt)
write.table(pleiotropy,file = "pleiotropy.txt",sep = "\t",quote = F)
heterogeneity=mr_heterogeneity(harm_rt)
write.table(heterogeneity,file = "heterogeneity.txt",sep = "\t",quote = F)
p1 <- mr_scatter_plot(mr_result, harm_rt)
ggsave(p1[[1]], file="scatter.pdf", width=8, height=8)
singlesnp_res<- mr_singlesnp(harm_rt)
singlesnpOR=generate_odds_ratios(singlesnp_res)
write.table(singlesnpOR,file="singlesnpOR.txt",row.names = F,sep = "\t",quote = F)
p2 <- mr_forest_plot(singlesnp_res)
ggsave(p2[[1]], file="forest.pdf", width=8, height=8)
sen_res<- mr_leaveoneout(harm_rt)
p3 <- mr_leaveoneout_plot(sen_res)
ggsave(p3[[1]], file="sensitivity-analysis.pdf", width=8, height=8)
res_single <- mr_singlesnp(harm_rt)
p4 <- mr_funnel_plot(singlesnp_res)
ggsave(p4[[1]], file="funnelplot.pdf", width=8, height=8)
presso=run_mr_presso(harm_rt,NbDistribution = 1000)
capture.output(presso,file = "presso.txt")

#多基因流程
foreach(i=genelist1, .errorhandling = "pass") %do%{
  expo_rt<- read_exposure_data(
    filename = paste0(i,"_expo_rt.txt"),
    sep = "\t",
    snp_col = "rsids2",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    eaf_col = "eaf",
    pval_col = "p",
    samplesize_col = "n")
  

  outc_rt <- read_outcome_data(
    snps = expo_rt$SNP,
    filename = "myoutcome.gz",
    sep = "\t",
    snp_col = "rsids",
    beta_col = "beta",
    se_col = "sebeta",
    effect_allele_col = "alt",
    other_allele_col = "ref",
    eaf_col = "af_alt",
    pval_col = "pval")
  
harm_rt <- harmonise_data(
  exposure_dat =  expo_rt, 
  outcome_dat = outc_rt,action=2)
harm_rt$R2 <- (2 * (harm_rt$beta.exposure^2) * harm_rt$eaf.exposure * (1 - harm_rt$eaf.exposure)) /
  (2 * (harm_rt$beta.exposure^2) * harm_rt$eaf.exposure * (1 - harm_rt$eaf.exposure) +
     2 * harm_rt$samplesize.exposure*harm_rt$eaf.exposure * (1 - harm_rt$eaf.exposure) * harm_rt$se.exposure^2)
harm_rt$f <- harm_rt$R2 * (harm_rt$samplesize.exposure - 2) / (1 - harm_rt$R2)
harm_rt$meanf<- mean( harm_rt$f)
harm_rt<-harm_rt[harm_rt$f>10,]

mr_result<- mr(harm_rt)
result_or=generate_odds_ratios(mr_result) 
dir.create(i) 
filename=i
write.table(harm_rt, file =paste0(filename,"/harmonise.txt"),row.names = F,sep = "\t",quote = F)
write.table(result_or[,5:ncol(result_or)],file =paste0(filename,"/OR.txt"),row.names = F,sep = "\t",quote = F)
pleiotropy=mr_pleiotropy_test(harm_rt)
write.table(pleiotropy,file = paste0(filename,"/pleiotropy.txt"),sep = "\t",quote = F)
heterogeneity=mr_heterogeneity(harm_rt)
write.table(heterogeneity,file = paste0(filename,"/heterogeneity.txt"),sep = "\t",quote = F)
p1 <- mr_scatter_plot(mr_result, harm_rt)
ggsave(p1[[1]], file=paste0(filename,"/scatter.pdf"), width=8, height=8)
singlesnp_res<- mr_singlesnp(harm_rt)
singlesnpOR=generate_odds_ratios(singlesnp_res)
write.table(singlesnpOR,file=paste0(filename,"/singlesnpOR.txt"),row.names = F,sep = "\t",quote = F)
p2 <- mr_forest_plot(singlesnp_res)
ggsave(p2[[1]], file=paste0(filename,"/forest.pdf"), width=8, height=8)
sen_res<- mr_leaveoneout(harm_rt)
p3 <- mr_leaveoneout_plot(sen_res)
ggsave(p3[[1]], file=paste0(filename,"/sensitivity-analysis.pdf"), width=8, height=8)
res_single <- mr_singlesnp(harm_rt)
p4 <- mr_funnel_plot(singlesnp_res)
ggsave(p4[[1]], file=paste0(filename,"/funnelplot.pdf"), width=8, height=8)
presso=run_mr_presso(harm_rt,NbDistribution = 1000)
capture.output(presso,file = paste0(filename,"/presso.txt"))
}



####共定位分析####

#装包
if(!require("remotes"))
  install.packages("remotes")
install.packages("dplyr")
library(remotes)
install_github("chr1swallace/coloc",build_vignettes=TRUE)
library("coloc")
library(dplyr)

#导入表型1数据
gwas <- read.table(file="E:/path_to_GWAS/GWAS.txt", header=T);
#导入表型2数据
eqtl <- read.table(file="E:/path_to_eqtl/eQTL.txt", header=T);

#数据合并
input <- merge(eqtl, gwas, by="rs_id", all=FALSE, suffixes=c("_eqtl","_gwas"))
head(input)

#共定位分析
result <- coloc.abf(dataset1=list(pvalues=input$pval_nominal_gwas, type="cc", s=0.33, N=50000), dataset2=list(pvalues=input$pval_nominal_eqtl, type="quant", N=10000), MAF=input$maf)

#筛选共定位位点
library(dplyr)
need_result=result$results %>% filter(SNP.PP.H4 > 0.95)#一般0.95





























































