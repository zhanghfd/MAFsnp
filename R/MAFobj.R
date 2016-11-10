MAFobj <-
function(d){
  count = d$count;
  S= ncol(count)/2;
  X = count[,1:S];N = count[,-c(1:S)];
  Xsum = rowSums(X);
  d$X = X;
  d$N = N;
  d$Xsum = Xsum;
  rm(count,X,N,Xsum);gc();
  return (d);
}
