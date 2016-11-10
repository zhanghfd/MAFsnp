SNP.calling <-
function(d,FDR=0.01){  
  qval = d$q.value;
  pval = d$p.value;
  Xsum = d$Xsum;
  M = length(Xsum);
  position = d$position[Xsum>0];
  
  called.snps = function(qval,index)
  {
    ind = which(qval<FDR);
    snp = index[ind];
    cbind(snp,pval[ind],qval[ind]);
  }
  SNPs = called.snps(qval,position);	
  colnames(SNPs)=c('position','p.value','q.value');
  d$snps = SNPs;
  rm(pval,qval,SNPs,position);gc();
  return (d);
}
