
chr = as.numeric(commandArgs(TRUE))
tissue = 'Whole Blood'


library(data.table)
.libPaths("~/R/rlib-3.4_new")
library('MendelianRandomization')



flip_dosage = function(x){
  if(x == 0){
    ans = 2
  }else if(x == 1){
    ans = 1
  }else if(x == 2){
    ans = 0
  }
}

#gene info
gene_info = as.data.frame(fread('/data/xxx/anno/gencode/37/gencode.v32.GRCh37.txt'))

#GTEx tissue sample info
sample_info = as.data.frame(fread('/data/xxx/data/gtex/exp/v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt'))
sample_info = sample_info[sample_info$SMTSD == tissue,]
#sample_info$SAMPID = gsub(pattern = '-', replacement = '.', x = sample_info$SAMPID)

#load GTEx APA
apa_all = as.data.frame(fread('/data/xxx/projects/meth_apa/apa_gtex/PDUI.txt'))
gene_list = apa_all$ensembl_gene
apa_all = as.data.frame(t(apa_all[,-c(1,2,3)]),stringsAsFactors = F)
colnames(apa_all) = gene_list
apa_all$sampleid = rownames(apa_all)
apa_all = apa_all[,c(ncol(apa_all), seq(1,(ncol(apa_all)-1)))]
apa_all = apa_all[which(apa_all$sampleid %in% sample_info$SAMPID),]
n_non_na = apply(apa_all, MARGIN = 2, function(x) length(which(!is.na(x))))
apa_all = apa_all[,which(n_non_na>=100)]   #filter n_APA<100
apa_var = apply(apa_all[,-1], MARGIN = 2, function(x) var(x,na.rm = T))
apa_all = apa_all[,c(1,which(apa_var!=0)+1)]   #filter var(APA)=0
apa_all$sampleid = sapply(apa_all$sampleid, function(x) paste0('GTEX-',strsplit(x,"[-]")[[1]][2]))  #sample id
apa_all = apa_all[!duplicated(apa_all$sampleid),]
geneid_non_na = which(colnames(apa_all) %in% gene_info$geneid)  #filter na gene name
apa_all = apa_all[,c(1,geneid_non_na)]


#load methylation annotation file
meth_anno = as.data.frame(fread('/data/xxx/anno/meth/HumanMethylation450_15017482_v1-2.csv',fill=TRUE))
meth_anno = meth_anno[,c('IlmnID','UCSC_RefGene_Name')]
colnames(meth_anno) = c('cpg','genename')

#load snp-meth association
load('/data/xxx/projects/meth_apa/snp_meth/raw/cis-cosmopairs_combined_151216.RData')
# > cosmo$
# cosmo$pair           cosmo$se.disco       cosmo$p.xe.disco
# cosmo$snp            cosmo$p.disco        cosmo$beta.xe.comb
# cosmo$snp.chr        cosmo$beta.repl      cosmo$se.xe.comb
# cosmo$snp.pos        cosmo$se.repl        cosmo$p.xe.comb
# cosmo$A1             cosmo$p.repl         cosmo$repli.xe
# cosmo$A2             cosmo$beta.comb      cosmo$pop.disco
# cosmo$eaf            cosmo$se.comb        cosmo$z
# cosmo$cpg            cosmo$p.comb         cosmo$tmpid
# cosmo$cpg.chr        cosmo$repli.pop      cosmo$tmpid2
# cosmo$cpg.pos        cosmo$beta.xe.disco
# cosmo$beta.disco     cosmo$se.xe.disco
snp_meth_pos = which(cosmo$cpg.chr == chr)
snp_meth = data.frame(rsid = cosmo$snp[snp_meth_pos],
                      cpg = cosmo$cpg[snp_meth_pos],
                      a1_snp_cpg = cosmo$A1[snp_meth_pos],
                      a2_snp_cpg = cosmo$A2[snp_meth_pos],
                      beta_snp_cpg = cosmo$beta.comb[snp_meth_pos],
                      se_snp_cpg = cosmo$se.comb[snp_meth_pos],
                      stringsAsFactors = F)
snp_meth$rsid = as.character(snp_meth$rsid)
snp_meth$cpg = as.character(snp_meth$cpg)

#load independent snps (LD r2=0.1)
snp_pruned = read.table('/data/xxx/data/gtex/geno/v8/prune0.1/v8_prune0.1_ea.bim',stringsAsFactors = F)

#need 3 snps for each cpg site as IV to perform MR
snp_meth = snp_meth[which(snp_meth$rsid %in% snp_pruned$V2),]
check_n_of_snp_per_cpg = as.data.frame(table(snp_meth$cpg),stringsAsFactors=F)
cpg_keep = check_n_of_snp_per_cpg[check_n_of_snp_per_cpg$Freq>=3,1]

#generate CpG-Gene df
info = data.frame(cpg = NA, genename = NA)

for(i in 1:length(cpg_keep)){
  gene_name_tmp = meth_anno[which(meth_anno$cpg == cpg_keep[i]),'genename']
  gene_name_array = sapply(gene_name_tmp, function(x) unique(strsplit(x,"[;]")[[1]]))
  if(length(gene_name_array[[1]])>=1){
    pos_start = nrow(info)+1
    pos_stop = nrow(info)+length(gene_name_array)
    info[pos_start:pos_stop,'genename'] = gene_name_array
    info[pos_start:pos_stop,'cpg'] = cpg_keep[i]
  }
}

info = merge(info, gene_info[,c('genename','geneid')], by='genename')
info = info[which(info$geneid %in% colnames(apa_all)),]

#perform MR for each CpG-Gene pair
tissue = gsub(pattern = ' ', replacement = '_', x = tissue)

for(i in 1:nrow(info)){
  cpg = info[i,'cpg']
  geneid = info[i,'geneid']
  
  #extract snps and convert to dosage
  snp_tmp = snp_meth[which(snp_meth$cpg == cpg),]

  write.table(snp_tmp,paste0('/data/xxx/projects/meth_apa/tmp/',tissue,'_',cpg,'_',geneid,'.genelist'), quote = F, sep='\t', row.names = F)
  
  cmd = paste0('plink --bfile /data/xxx/data/gtex/geno/v8/prune0.1/v8_prune0.1_ea --extract /data/xxx/projects/meth_apa/tmp/',tissue,'_',cpg,'_',geneid,'.genelist --recode A --out /data/xxx/projects/meth_apa/tmp/',tissue,'_',cpg,'_',geneid)
  system(cmd)
  
  #load genotype in dosage
  dosage = read.table(paste0('/data/xxx/projects/meth_apa/tmp/',tissue,'_',cpg,'_',geneid,'.raw'), header = T, stringsAsFactors = F)
  dosage = dosage[,-c(1,3,4,5,6)]
  
  #df for MR
  df = data.frame(rsid = sapply(colnames(dosage[-1]), function(x) strsplit(x,"[_]")[[1]][1]),
                  a1_apa = sapply(colnames(dosage[-1]), function(x) strsplit(x,"[_]")[[1]][2]),
                  snp_apa_beta = NA,
                  snp_apa_se = NA)
  
  tmp = merge(apa_all[,c(1,which(colnames(apa_all)==geneid))],dosage,by=1)
  
  #estimate SNP-APA association
  for(j in 1:nrow(df)){
    snp = df[j,'rsid']
    ans = summary(lm(tmp[,geneid]~tmp[,j+2]))
    df[j,'snp_apa_beta'] = ans$coefficients[2,1]
    df[j,'snp_apa_se'] = ans$coefficients[2,2]
  }
  
  #MR
  df = merge(df, snp_tmp, by='rsid')
  df$snp_apa_beta = ifelse(df$a1_apa == df$a1_snp_cpg, df$snp_apa_beta, df$snp_apa_beta*-1) #flip allele
  
  #flip dosage allele
  colnames(dosage) = sapply(colnames(dosage), function(x) strsplit(x,"[_]")[[1]][1])
  dosage = dosage[,-1]
 
  for(flip_i in 1:nrow(df)){
    if(df$a1_apa[flip_i] != df$a1_snp_cpg[flip_i]){
      
      if(length(which(is.na(dosage[,flip_i])))!=0){
        dosage[which(is.na(dosage[,flip_i])),flip_i] = median(dosage[,flip_i], na.rm = T) 
      }
      dosage[,flip_i] = sapply(dosage[,flip_i], function(x) flip_dosage(x))
    }
  }
  
  snp_correlation = cor(dosage)
  
  
  
  
  #----MR-Weighted_median----
  set.seed(as.numeric(sub('ENSG','',geneid))+100)
  ans_median_weighted<-try(mr_median(mr_input(bx = df$beta_snp_cpg, 
                                          bxse = df$se_snp_cpg, 
                                          by = df$snp_apa_beta, 
                                          byse = df$snp_apa_se,
                                          correlation = snp_correlation),
                                 weighting = "weighted", iterations = 100))
  if('try-error' %in% class(ans_median_weighted)){  
    next
  }
  
  info[i,'n_IVs']<-nrow(df)
  info[i,'median_weighted_beta']<-ans_median_weighted@Estimate
  info[i,'median_weighted_se']<-ans_median_weighted@StdError
  info[i,'median_weighted_p']<-ans_median_weighted@Pvalue
  
  #if p<0.05, save the MR df
  if(ans_median_weighted@Pvalue<0.05){
    write.table(df, paste0('/data/xxx/projects/meth_apa/snp_meth_apa_mr/mr_df/',tissue,'_',cpg,'_',geneid,'.txt'),quote = F,sep='\t',row.names = F)
  }
  
}


write.table(info,paste0('/data/xxx/projects/meth_apa/snp_meth_apa_mr/tmp/',tissue,'_',chr,'.txt'),quote = F,sep = '\t',row.names = F)



