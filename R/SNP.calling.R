SNP.calling <-
function(d,FDR=0.01){  
  p.fdr = d$p.fdr;
  p.val = d$p.value;
  Xsum = d$Xsum;
  M = length(Xsum);
  position = d$position[Xsum>0];
  
  called.snps = function(p.fdr,index)
  {
    ind = which(p.fdr<FDR);
    cbind(index[ind],p.fdr[ind],p.fdr[ind]);
  }
  SNPs = called.snps(p.fdr,position);	
  colnames(SNPs)=c('position','p.value','p.fdr');
  d$snps = SNPs;
  rm(p.fdr,SNPs,position);
  gc();
  return (d);
}
