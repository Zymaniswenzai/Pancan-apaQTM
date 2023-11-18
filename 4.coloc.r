library(coloc)
cancer=c('BRCA','CESC','COAD','ESCA',
         'HNSC','KIRC','KIRP','LIHC',
         'LAML','LGG','LUAD','LUSC','PAAD',
         'PRAD','STAD','THCA','UCEC')
dd=c("breast_cancer_de_cpg.txt",
     "cervical_cancer_de_cpg.txt",
     "colorectal_cancer_de_cpg.txt",
     "esophageal_carcinoma_de_cpg.txt",
     "head_and_neck_squamous-cell_carcinoma_de_cpg.txt",
     "clear_cell_renal_cell_carcinoma_de_cpg.txt",
     "kidney_renal_papillary_cell_carcinoma_de_cpg.txt",
     "hepatocellular_carcinoma_de_cpg.txt",
     "acute_myeloid_leukemia_de_cpg.txt",
     "glioma_de_cpg.txt",
     "lung_adenocarcinoma_de_cpg.txt",
     "lung_squamous_cell_carcinoma_de_cpg.txt",
     "pancreatic_cancer_de_cpg.txt",
     "prostate_cancer_de_cpg.txt",
     "stomach_cancer_de_cpg.txt",
     "thyroid_cancer_de_cpg.txt",
     "uterine_carcinosarcoma_de_cpg.txt")

wd=c('/home/lailab/disk/lyg/Kb100/apaMeth/' ,
     '/home/lailab/disk/lyg/Kb100/expMeth/',
     '/home/lailab/disk/lyg/Kb100/splMeth/',
     '/home/lailab/disk/zdd/APA/EWAS12/PEERedAPAsd/',
     '/home/lailab/disk/zdd/APA/EWAS12/PEERedEXPsd/',
     '/home/lailab/disk/zdd/APA/EWAS12/PEERedSplicesd/')
for (i in 1:length(cancer)){
  
  ewas <- read.table(file=dd[i], header=T,sep='\t')
  ewas=ewas[which(ewas$beta!=''),]
  
  ewas$varbeta=ewas$se^2
  ewas$CpG=rownames(ewas)
  
  sd=read.csv(paste0(cancer[i],'_sd.csv'))
  file1=paste0(cancer[i],'100kbCis.txt')
  data=read.table(file1,header=T,sep='\t')
  
  data$se=data$beta/data$statistic
  data$varbeta=data$se^2
  
  apa=sd$id
  std=data$gene
  ustd=unique(std)
  res=c()
  length(unique(data$gene))
  for(j in 1:length(unique(data$gene))){
    sdY=as.numeric(sd$sd[which(apa==ustd[j])])
    ss=which(std==ustd[j])
    snp=data$snps[ss]
    pvalue=data$pvalue[ss]
    beta=data$beta[ss]
    se=data$beta[ss]/data$statistic[ss]
    varbeta=se^2
    dd1=list(snp=snp,sdY=sdY,pvalues=pvalue,beta=beta,varbeta=varbeta,type="quant")
    
    ee=ewas[which(!duplicated(ewas$CpG)),]
    snp=ee$CpG
    pvalue=ee$p
    verbeta=ee$varbeta
    beta=ee$beta
    dd2=list(snp=snp,type="cc",varbeta=verbeta,beta=beta)
    sn=intersect(dd1$snp,dd2$snp)
    if(length(sn)>0){
      re=try( coloc.abf(dataset1 = dd1,dataset2 = dd2),silent =TRUE)
      if(length(re$results)>10){
        re1=re$results
        re1$PP3=re$summary[[5]]
        re1$PP4=re$summary[[6]]
      }
      if(length(re1)>10){
        re1=cbind(re1,ustd[j])
        res=rbind(res,re1)
      }
    }
  }
  
  write.table(res,paste0(cancer[i],'_coloc100kb.txt'),sep='\t',row.names = F,quote=F)
} 
