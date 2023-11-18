library(mediation)

aqtm=read.csv('BRCATransSig.csv',row.names = 1)
eqtm=read.table('BRCA_CisSig.txt',sep = '\t',header = T)
aparegu=read.table('APAregu.txt',sep = '\t',header = T)
dat=read.table('BRCA_meth_APA.txt',sep = '\t',header = T,check.names = F)
exp=read.table('BRCA_APAregulator_expression.txt',header = T,check.names = F)

aqtm=aqtm[,1:2]
eqtm=eqtm[,1:2]
eqtm=merge(eqtm,aparegu,by='gene')
aqtm=merge(eqtm,aqtm,by='snps')

res=c()
rtmp=data.frame(cpg=NA,APA=NA,gene=NA,ACME_Estimate=NA,
                ACME_pval=NA,ACME_CI_lower=NA,ACME_CI_Upper=NA,
                ADE_Estimate=NA,ADE_pval=NA,ADE_CI_lowerr=NA,
                ADE_CI_Upperr=NA,   total_effect_Estimater=NA,
                total_effect_pvalr=NA,TE_CI_lower=NA,
                TE_CI_Upper=NA, Prop.Mediated_Estimate=NA,
                Prop.Mediated_pval=NA,Prop_CI_lower=NA, 
                Prop_CI_Upper=NA)
for(i in 1:dim(aqtm)[1]){
  atmp=aqtm[i,3]
  ctmp=aqtm[i,1]
  tmp=dat[,which(colnames(dat)==atmp|colnames(dat)==ctmp)]
  for(j in 1:ncol(exp)){
    tmp1=cbind(tmp,exp[j])
    names(tmp1)=c('Predictor','Outcome','Mediator')
    model_mediator <- lm(Mediator ~ Predictor, data = tmp1)
    model_outcome <- lm(Outcome ~ Predictor + Mediator, data = tmp1)
    set.seed(123.234)
    mediation_result <- mediate(model_mediator,model_outcome,treat ='Predictor',mediator = 'Mediator',bootstrap = 1000)
    rr1=summary(mediation_result)
    rtmp[j,1]=names(tmp)[1]
    rtmp[j,2]=names(tmp)[2]
    rtmp[j,3]=names(exp)[j]
    rtmp[j,4]=rr1$d0
    rtmp[j,5]=rr1$d0.p
    rtmp[j,6]=rr1$d0.ci[1]
    rtmp[j,7]=rr1$d0.ci[2]
    rtmp[j,8]=rr1$z0
    rtmp[j,9]=rr1$z0.p
    rtmp[j,10]=rr1$z0.ci[1]
    rtmp[j,11]=rr1$z0.ci[2]
    rtmp[j,12]=rr1$tau.coef
    rtmp[j,13]=rr1$tau.p
    rtmp[j,14]=rr1$tau.ci[1]
    rtmp[j,15]=rr1$tau.ci[2]
    rtmp[j,16]=rr1$n.avg
    rtmp[j,17]=rr1$n.avg.p
    rtmp[j,18]=rr1$n.avg.ci[1]
    rtmp[j,19]=rr1$n.avg.ci[2]
  }
  rtmp1=rtmp[which(rtmp$ACME_pval<0.05),]
  rtmp1=rtmp1[which(rtmp1$ACME_Estimate!=''),]

  if(dim(rtmp1)[1]>0){
    res=rbind(res,rtmp1)
  }else{i=i+1}
}

write.table(res,'BRCA_mediation_result.txt',sep='\t',quote=F,row.names = F)
