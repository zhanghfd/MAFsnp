eLRT <-
function(d){
  position = d$position;
  eps = 1e-10;
  X= d$X;
  N = d$N;
  if(length(d$p)>0){
    p.true = d$p;
    p.true[p.true<eps] = eps;
  }else{
    p.true = rep(nrow(X),0.01);
  }
  Xsum = d$Xsum;
  M = length(Xsum);
  ind = which(Xsum>0);
  X = X[ind,];N = N[ind,];p.true = p.true[ind];
  dbinom.vec = Vectorize(dbinom,vectorize.args = c("x","size","prob"),SIMPLIFY = FALSE);
  optim.func = function(t1,t2,t3,p0,try2=FALSE){
    f1= function(p){
      ps = c(p^2,2*p*(1-p),(1-p)^2);	
      ps[1]* t1+ ps[2]* t2+ps[3]* t3;
    }
    fn1 = function(p){ 		
      -sum(log(f1(p)));		
    }
    gn1 = function(p){ 
      ps = c(2*p,2*(1-p),-2*(1-p));
      g1 = ps[1]* t1+ ps[2]* t2+ps[3]* t3;
      -sum(g1/f1(p));		
    }
    
    lower = 1e-20;upper = 0.5-1e-4;
    solve = optim(p0,fn1,gn1,method='L-BFGS-B',lower=lower,upper=upper,control=list(maxit = 1e4));
    convergence =solve$convergence;t=1;
    if(try2){
      while(convergence!=0 & t<40){
        th = runif(1,lower,upper);
        solve = optim(th,fn1,gn1,method='L-BFGS-B',lower=lower,upper=upper,control=list(maxit = 1e4));
        convergence = solve$convergence;
        t = t+1;
      }
    }
    return (c(as.numeric(unlist(solve)),t));		
  }
  
  optim.vec = Vectorize(optim.func,vectorize.args = c("t1","t2","t3","p0"),SIMPLIFY = TRUE);
  M.pos = sum(Xsum>0);
  m = 5e3; n = ifelse(M.pos%%m ==0,M.pos%/%m,M.pos%/%m+1);
  res = matrix(NA,nrow = M.pos,ncol = 7);
  colnames(res)=c('par','value','counts.function','counts.gradient','convergence','message','t');
  e.hat = fn0s= rep(NA,M.pos);
  for (i in 1:n){
    s = (i-1)*m;
    if(i <n){
      sub = 1:m + s;
    }else{
      sub = (s+1):M.pos;
    }
    Xsub = X[sub,];
    Nsub = N[sub,];
    X.list = unclass(as.data.frame(t(Xsub)));
    N.list = unclass(as.data.frame(t(Nsub)));
    e.hat.sub = sapply(X.list,mean)/sapply(N.list,mean);
    p0 = p.true[sub];   
    k = which(e.hat.sub>0.5);
    if(length(k)!=0){
      X.list[k]= mapply(function(x,n){n-x},X.list[k],N.list[k],SIMPLIFY=FALSE)
    }
    
    j = which(e.hat.sub>0.1);
    p0[j] = e.hat.sub[j];
    a = which(p0>=0.5);p0[a]=1-p0[a];    
    e.hat.sub[j] = 0.05;
    e.hat[sub] = e.hat.sub;
    T1 = dbinom.vec(x = X.list,size = N.list,prob = 1-e.hat.sub);
    T2 = dbinom.vec(x = X.list,size = N.list,prob = .5);	
    T3 = dbinom.vec(x = X.list,size = N.list,prob = e.hat.sub);
    fn0s.sub = -colSums(log(as.matrix(as.data.frame(T3)))); #vector.
    fn0s[sub] = fn0s.sub; 
    res1 = t(optim.vec(T1,T2,T3,p0));
    unconvergence = which(res1[,5]!=0); 
    if(length(unconvergence)>0){
      res2 = t(optim.vec(T1[unconvergence],T2[unconvergence],T3[unconvergence],p0[unconvergence],try2=TRUE));		
      res1[unconvergence,]=res2;	
    }
    res[sub,] = res1;
  }	
  psudo.lrt = 2*(fn0s-res[,'value']);	
  d$position.with.variation = d$position[d$Xsum>0];
  d$T = psudo.lrt;
  rm(psudo.lrt,T1,T2,T3,X.list,N.list);gc();
  return(d);
  
}
