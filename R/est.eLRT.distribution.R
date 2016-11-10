est.eLRT.distribution <-
function(d,cutoff = 15){
  flag = d$flag[d$Xsum>0];
  T = d$T;  
  Tnull = T[flag==0]; 
  chisq0 = rchisq(1e4,df = 1);  
  ind = which(Tnull<1e-8); l = length(ind);
  a = l/length(Tnull);  
  Tnull.positive = Tnull[-ind];
  cutoff = min(cutoff,max(Tnull.positive,na.rm=TRUE));
  chisq0 = chisq0[chisq0<cutoff];	
  Tnull.positive=Tnull.positive[Tnull.positive<cutoff];
  k = mean(Tnull.positive,na.rm=TRUE)/mean(chisq0);
  a.k = list(a.est=a,k.est = k);
  d$a.k = a.k;
  rm(flag,T,Tnull,chisq0,Tnull.positive);gc();
  return (d);
}
