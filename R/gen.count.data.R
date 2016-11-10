gen.count.data <-
function(ErrorProb,nSample,N,SNPtype,M=1e3,r=100){
  mean.by.var = 0.4;
  ps = c(SNPtype,rep(0,r));
  P.all =ERROR.all=FLAG.all= rep(NA,M*(r+1));
  GENOTYPE.all= matrix(NA,ncol=3,nrow = M*(r+1)); 
  COUNT.all = matrix(NA,ncol=nSample*2,nrow = M*(r+1));
  
  get.p = function(p,M){
    P = NULL;
    if(p==0){
      P = rep(0,M); 
    }else{
      if(p==1){
        P = runif(M,0.001,0.01); 
      }else{
        if(p==2){
          P = runif(M,0.01,0.05);
        }else{ 
          P = runif(M,0.05,0.1);        
        }
      }
    }
    return (P);
  }
  get.genotype = function(M,nSample,Ps){
    genotype = NULL;
    Ps = unclass(as.data.frame(t(Ps))); 
    vec.rmultinom = Vectorize(rmultinom,"prob");
    genotype = t(vec.rmultinom(1,nSample,Ps));
    colnames(genotype)= c('RR','Rr','rr');
    return (genotype);
  }
  get.coverage = function(M,nSample,N,is.pois,mean.by.var){
    R.num = NULL;
    R.num = rgenpois(nSample*M,N,mean.by.var=mean.by.var);
    R.num = matrix(R.num,ncol = nSample);
    return (R.num);
  }
  get.e = function(M,e){
    error = rnorm(2*M,e,sd = 1e-3);
    error = error[error>0][1:M];
    error;
  }
  get.X = function(genotype,coverage,error,nSample){
    M = dim(genotype)[1];
    X = NULL; 
    S = nSample;
    coverage = unclass(as.data.frame(coverage)); 
    genotype = unclass(as.data.frame(t(genotype)));
    genotype =t(sapply(genotype,function(g){c(rep(0,g[1]),rep(1,g[2]),rep(2,g[3]))}));
    
    probs = matrix(NA,nrow = M,ncol=S);
    
    for (s in 1:S){  
      g = genotype[,s];
      p = rep(.5,M);
      i0 = which(g==0);p[i0]=error[i0];
      i2 = which(g==2);p[i2]=1-error[i2];		
      probs[,s]=p;
    }	
    probs = unclass(as.data.frame(probs));
    vec.rbinom = Vectorize(rbinom,vectorize.args = c("size","prob"));
    X = mapply(vec.rbinom,size = coverage,prob = probs,MoreArgs = list(n=1));
    rm(coverage,genotype,probs);gc();	
    return (X); 
  }
  rgenpois = function(Num,N,mean.by.var){
    ld = 1-sqrt(mean.by.var);
    th = N * (1-ld);
    n.max = 500;
    tmp = log(1:n.max);
    tmp = c(0,cumsum(tmp));
    rs = th*exp((0:n.max-1)*log(th+(0:n.max)*ld) - th - (0:n.max)*ld - tmp);
    rs = rs / sum(rs);
    res = sample(0:n.max,Num,prob=rs,replace=TRUE);
    return(res);
  }  
  
  for (i in 1:(r+1)){
    start = (i-1)*M +1;end = i*M;
    ind = start:end;
    p= ps[i];
    P = get.p(p,M); 
    P.all[ind] = P;
    Ps = cbind((1-P)^2, 2*P*(1-P), P^2); 
    GENOTYPE = get.genotype(M,nSample,Ps);
    FLAG = ifelse(rowSums(GENOTYPE[,2:3])!=0,1,0); 
    FLAG.all[ind]=FLAG;
    COVERAGE = get.coverage(M,nSample,N,mean.by.var=mean.by.var)	
    ERROR = get.e(M,ErrorProb); 
    ERROR.all[ind]= ERROR;
    VARIANT= get.X(GENOTYPE,COVERAGE,ERROR,nSample);
    GENOTYPE.all[ind,] = GENOTYPE;  
    count.MAFsnp = cbind(VARIANT,COVERAGE);	
    COUNT.all[ind,] = count.MAFsnp;
  }
  position = 1:(M*(r+1));
  colnames(COUNT.all)=c(paste("X",1:nSample,sep=''),paste("N",1:nSample,sep=''));
  colnames(GENOTYPE.all)= c('RR','Rr','rr');
  dat = list(count = COUNT.all,p=P.all,genotype=GENOTYPE.all,error = ERROR.all,flag=FLAG.all,position=position);	
  rm(GENOTYPE.all,ERROR.all,P.all,FLAG.all,COUNT.all);
  gc();
  return(dat);
}
