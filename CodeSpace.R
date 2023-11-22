# Generate PQX random variables by MC method.
rPQX<-function(alpha, theta, n){
  u<-runif(n,0,1)
  lambda<-sapply(u, function(x){ifelse(x<=alpha/(alpha+1), rexp(1, rate=theta), rgamma(1, shape=3, rate=theta))})
  res<-sapply(lambda,  function(x) rpois(1,lambda = x))
  return(res)
}

# My simulated data with parameter alpha = 1 and theta=0.5
X<-rPQX(24.89,0.264, 50000)
plot_ly()%>%add_histogram(x=~X, type = "histogram", histnorm = "probability density", nbinsx=150, marker=list(color="skyblue", line=list(color="black", width=1.5)))%>%
  add_trace(x=density(X)$x, y=density(X)$y, type="scatter", mode="lines", line=list(color="blue", width=1.5))%>%
  layout(xaxis=list(range=c(-0.5,20.5)))



# Bootstrap estimation
BootX<-sapply(1:10000, function(i){
  samp<-sample(X, size = length(X), replace = TRUE)
  m1<-mean(samp)
  m2<-samp%*%samp/length(samp)
  a<- (-7*m1^2+sqrt(25*m1^4+12*m1^3-12*m1^2*m2) -3*(m1-m2))/(2*m1^2+m1-m2)
  b<- (a+3)/(m1*(1+a))
  return(c(a, b))
})
mean(BootX[1,])
mean(BootX[2,])



# Time Series
alpha<-0.35
theta<-0.15
p<-0.3

X<-rep(1, 1000)
for (i in 2:1000) {
  X[i]<-  sum(rbinom(X[i-1], 1, p))+rPQX(alpha, theta, 1)
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


#Real data
eqarchive<-read.csv("eqarchive-en.csv")
eqarchive<-eqarchive%>%mutate(YM=paste0(year(date), "-", month(date), "-1"))
data<-eqarchive%>%mutate(YM=as.POSIXct(YM))%>%filter(magnitude>=4)%>%group_by(YM)%>%summarise(n=n())
YW(data=data$n)
MLE(data=data$n, p=c(0.04,0.02,0.3))

EqCount_pred<-rep(data$n[1], length(data$n))
for (i in 2:length(data$n)) {
  EqCount_pred[i]<-  sum(rbinom(EqCount_pred[i-1], 1, 0.3666483))+rPQX(24.8872548, 0.2637256,1)
}
data.frame(Time=data$YM, EqCount=data$n, EqCount_pred)%>%gather(key = "Type",value = "value", EqCount:EqCount_pred)%>%
  plot_ly(x=~Time, y=~value, color = ~Type, type="scatter", mode="lines",
          line=list(color=~c(EqCount="black", EqCount_pred="red")[Type]))%>%
  layout(paper_bgcolor="lightgrey", plot_bgcolor="lightgrey")








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
res<-Sim_INARPQX$new(size=1000, nboot=100, plotMSE=FALSE)
res<-Sim_INARPQX$new(size=2000, nboot=100, plotMSE=TRUE)
res<-Sim_INARPQX$new(size=2000, nboot=100, plotMSE=TRUE, initialp=c(0.04,0.02,0.3))

res$estimates

YW_Boot(data=data, 1000000)


