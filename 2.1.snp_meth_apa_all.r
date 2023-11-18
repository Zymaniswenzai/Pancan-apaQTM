
args = as.numeric(commandArgs(TRUE))  #677
tissue = 'Whole Blood'


library(data.table)
.libPaths("~/R/rlib-3.4_new")
library('MendelianRandomization')
# library(foreach)
# library(doParallel)

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


#load independent snps (LD r2=0.1)
snp_pruned = read.table('/data/xxx/data/gtex/geno/v8/prune0.1/v8_prune0.1_ea.bim',stringsAsFactors = F)

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

#all cpgs
cpg_list_all = as.character(unique(cosmo$cpg))
i_start = (args-1)*100+1
i_end = min(args*100,length(cpg_list_all))

cpg_list = cpg_list_all[i_start:i_end]

cpg_i = 3
for(cpg_i in 1:length(cpg_list)){
  print(cpg_i)
  



# registerDoParallel(cores = detectCores())
# 
# trials <- length(gene_list)
# 
# system.time({
#   result <- foreach(i = 1:trials, .combine=rbind) %dopar% {
#     gene_id = gene_list[i]
#     unlist(asso(gene_id))
#   }
# })
# 
# stopImplicitCluster()




cpg = cpg_list[cpg_i]

snp_meth_pos = which(cosmo$cpg == cpg)


snp_meth = data.frame(rsid = cosmo$snp[snp_meth_pos],
                      cpg = cosmo$cpg[snp_meth_pos],
                      a1_snp_cpg = cosmo$A1[snp_meth_pos],
                      a2_snp_cpg = cosmo$A2[snp_meth_pos],
                      beta_snp_cpg = cosmo$beta.comb[snp_meth_pos],
                      se_snp_cpg = cosmo$se.comb[snp_meth_pos],
                      stringsAsFactors = F)
#rm(cosmo)

snp_meth$rsid = as.character(snp_meth$rsid)
snp_meth$cpg = as.character(snp_meth$cpg)

snp_meth = snp_meth[which(snp_meth$rsid %in% snp_pruned$V2),]

print(nrow(snp_meth))
if(nrow(snp_meth)<3){next}

#generate CpG-Gene df
info = data.frame(cpg = cpg, geneid = colnames(apa_all)[-1], stringsAsFactors = F)


#extract snps and convert to dosage
tissue = gsub(pattern = ' ', replacement = '_', x = tissue)
write.table(snp_meth,paste0('/data/xxx/projects/meth_apa/tmp/',tissue,'_',cpg,'.genelist'), quote = F, sep='\t', row.names = F)
cmd = paste0('plink --bfile /data/xxx/data/gtex/geno/v8/prune0.1/v8_prune0.1_ea --extract /data/xxx/projects/meth_apa/tmp/',tissue,'_',cpg,'.genelist --recode A --out /data/xxx/projects/meth_apa/tmp/',tissue,'_',cpg)
system(cmd)


#load genotype in dosage
dosage = read.table(paste0('/data/xxx/projects/meth_apa/tmp/',tissue,'_',cpg,'.raw'), header = T, stringsAsFactors = F)
dosage = dosage[,-c(1,3,4,5,6)]

#rm tmp file
cmd = paste0('rm /data/xxx/projects/meth_apa/tmp/',tissue,'_',cpg,'*')
system(cmd)

sample_list = intersect(dosage$IID, apa_all$sampleid)

dosage = dosage[sapply(sample_list, function(x) which(dosage$IID == x)),]
apa = apa_all[sapply(sample_list, function(x) which(apa_all$sampleid == x)),]
dosage = dosage[,-1]

df_raw = data.frame(rsid = sapply(colnames(dosage), function(x) strsplit(x,"[_]")[[1]][1]),
                a1_apa = sapply(colnames(dosage), function(x) strsplit(x,"[_]")[[1]][2]))

for(i in 1:nrow(info)){
  #print(i)
  geneid = info[i,'geneid']
  df = df_raw
  
  for(j in 1:ncol(dosage)){
    
    ans = summary(lm(apa[,geneid]~dosage[,j]))
    if(nrow(ans$coefficients)<2){
      df[j,'snp_apa_beta'] = 0
      df[j,'snp_apa_se'] = 99
    }else if(ans$coefficients[2,2] == 0){
      df[j,'snp_apa_beta'] = 0
      df[j,'snp_apa_se'] = 99
    }else{
      df[j,'snp_apa_beta'] = ans$coefficients[2,1]
      df[j,'snp_apa_se'] = ans$coefficients[2,2]
    }

    
  }
  
  df = merge(df, snp_meth, by = 'rsid')
  
  #MR
  df$snp_apa_beta = ifelse(df$a1_apa == df$a1_snp_cpg, df$snp_apa_beta, df$snp_apa_beta*-1) #flip allele
  
  #flip dosage allele
  colnames(dosage) = sapply(colnames(dosage), function(x) strsplit(x,"[_]")[[1]][1])
  
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
  
}



write.table(info,paste0('/data/xxx/projects/meth_apa/snp_meth_apa_mr/tmp/',tissue,'_',cpg,'.txt'),quote = F,sep = '\t',row.names = F)




}



















