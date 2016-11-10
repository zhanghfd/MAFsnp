pvalues <-
function(d){
  a.k = d$a.k;
  a = a.k$a.est;k = a.k$k.est;
  T = d$T;
  p.new = q.new = rep(NA,length(T));
  ind = which(T<1e-10);
  p.new[-ind] = 1- a - (1-a)*pchisq(T[-ind]/k,df=1);
  p.new[ind]=1;
  q.new = p.adjust(p.new,'fdr');
  
  d$p.value = p.new;
  d$q.value = q.new;
  rm(p.new,q.new,T);gc();
  return (d);
}
