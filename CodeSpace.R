# Generate PQX random variables by MC method.
rPQX<-function(alpha, theta, n){
  u<-runif(n,0,1)
  lambda<-sapply(u, function(x){ifelse(x<=alpha/(alpha+1), rexp(1, rate=theta), rgamma(1, shape=3, rate=theta))})
  res<-sapply(lambda,  function(x) rpois(1,lambda = x))
  return(res)
}

# YW
YW<-function(data){
  m<-mean(data)
  pb<-data.frame(data, L=lag(data), M=m)%>%na.omit()%>%mutate(numerator=(data-M)*(L-M))
  pb<-sum(pb$numerator)/((data-m)%*%(data-m))
  
  DI<-var(data)/m
  
  alpha<- (-3*DI+4*m+pb*(-3*DI-4*m+3)-sqrt(m*(pb-1)*(12*DI-13*m+pb*(12*DI+13*m-12)-12))+3)/(DI-m+pb*(DI+m-1)-1)
  theta<- (alpha+3)/(m*(1+alpha)*(1-pb))
  return(c(alpha, theta, pb))
}

# MLE
MLE<-function(data, p){
  l<-function(p, x){
    alpha<-p[1]
    theta<-p[2]
    pb<-p[3]
    res<-0
    
    for(t in 2:length(x)){
      total<-0
      for(i in 1:min(x[t], x[t-1])){
        part1<-dbinom(i, size=x[t-1], prob=pb)
        part2<-2*alpha*theta*(theta+1)^2 + theta^3*(x[t]-i+1)*(x[t]-i+2)
        part3<-2*(alpha+1)*(theta+1)^(x[t]-i+3)
        total<-total+part1*ifelse(alpha>0 & theta>0, 1,0)*part2/part3*ifelse( pb>0 & pb<1, 1,0)
      }
      res<-res+log(total)
    }
    return(-res)
  }
  
  nlm(f=l, p=p, x=data)
}


Sim_INARPQX<-setRefClass(
  "BootClass",
  fields = list(size="numeric",nboot="numeric",BootData="matrix",estimates="matrix", MSE="list",plotMSE="logical", initialp="ANY"),
  methods=list(
    initialize=function(nboot, plotMSE,size=NULL, initialp=NULL){
      if(!is.null(size)){.self$size<-size}else{stop("Must have the size of the sample specified.")}
      
      if(plotMSE==TRUE){
        if(!is.null(initialp) & class(initialp)=="numeric"){.self$initialp<-initialp}else{.self$initialp<-NULL}
        .self$MSE<-.self$SimulateMSE(.self$size, nboot, .self$initialp)}
      else{
        .self$BootData<-.self$TakeSample(.self$size, nboot)
        .self$estimates<-.self$YWcalculator(BootData)}
    }, 
    
    TakeSample=function(size, nboot){
      alpha<-0.35
      theta<-0.15
      p<-0.3
      
      NewData<-sapply(1:nboot, function(n){
        X<-rep(1, size)
        for (i in 2:size) {
          X[i]<-  sum(rbinom(X[i-1], 1, p))+rPQX(alpha, theta, 1)
        }
        return(X)
      })
      return(NewData)
    },
    
    YWcalculator=function(NewData){
      res<-apply(NewData,2, function(x) YW(data=x))
      cat("YW Parameter estimates: ", rowMeans(res), "\n")
      return(res)
    },
    
    MLEcalculator=function(NewData, initialp){
      res<-apply(NewData,2, function(x) {y<-MLE(data=x, p=initialp); return(y$estimate)})
      cat("MLE Parameter estimates: ", rowMeans(res), "\n")
      return(res)
    },
    
    SimulateMSE=function(size, nboot, initialp){
      cat("Calculating MSE and generating plot ---------------","\n")
      n<-seq(200,size,length.out=ceiling(size/200))%>%round()
      
      ests<-lapply(n, function(s){
        samp<-.self$TakeSample(s, nboot)
        
        if(is.null(initialp)==TRUE){
          df<-.self$YWcalculator(samp)%>%t()%>%data.frame(n=s)}
        else{
          df<-.self$MLEcalculator(samp,initialp)%>%t()%>%data.frame(n=s)}
        colnames(df)<-c("alpha_est", "theta_est", "p", "n")
        
        cat("Process completed --- ", round(s/max(n)*100,4), "%", "\n")
        return(df)
      })
      
      mse<-map(ests, function(df){
        df%>%transmute(alphaMSE=(alpha_est-0.35)^2, theta_MSE=(theta_est-0.15)^2, n)%>%colMeans()
      })
      
      do.call(rbind, mse)%>%data.frame()%>%
        gather("param", "value", 1:2)%>%
        plot_ly(x=~n, y=~value, color = ~param, type = "scatter", mode="lines")%>%print()
      return(mse)
    }
  )
)
#res<-Sim_INARPQX$new(size=1000, nboot=100, plotMSE=FALSE)
#res<-Sim_INARPQX$new(size=2000, nboot=100, plotMSE=TRUE)
#res<-Sim_INARPQX$new(size=2000, nboot=100, plotMSE=TRUE, initialp=c(0.04,0.02,0.3))
