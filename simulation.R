


#-----------------------------------------------------------------------------------------------------#
#   Simulation Study of 95% HDI Coverage Probability of F1 Scores Using Latent Class Analysis         #
#-----------------------------------------------------------------------------------------------------#



### Notes.
# The following is an example simulation script, with the following values fixed:
# 1) Prevalence is 0.1
# 2) Sample size is set at 50
# 3) All hyperparameters of the prior distribution are Jeffreys prior 

library(HDInterval);library(caret) ;library(BayesLCA)



## sample size
n<-50

## prevalence
pre<-0.1 

## Hyperparameter of Se
alpha_se<-1/2 
beta_se<-1/2

## Hyperparameter of PPV
alpha_ppv<-1/2  ## Hyperparameter
beta_ppv<-1/2   ## Prevalence


cov<-NULL
ppvdata<-NULL
sedata<-NULL
preppv<-NULL
prese<-NULL
af1<-NULL
bf1<-NULL
ppv<-NULL
se<-NULL
bf2<-NULL

for(i in 1:10000){
  ds<-rbinom(n,1,pre)
  
  
  ds1<-subset(ds,ds==1)
  ds0<-subset(ds,ds==0)
  
  
  
  A<-rbinom(length(ds1),1,0.9)
  B<-rbinom(length(ds1),1, 0.8)
  C<-rbinom(length(ds1),1, 0.7)
  ds<-data.frame(A,B,C)
  ds$pred<-apply(ds, 1, sum)
  ds$pred <- ifelse(ds$pred >= 2, 1,0)
  
  
  
  ds1<-data.frame(ds,ds1)
  
  A<-rbinom(length(ds0),1,0.3)
  B<-rbinom(length(ds0),1,0.2)
  C<-rbinom(length(ds0),1,0.1)
  ds<-data.frame(A,B,C)
  ds$pred<-apply(ds, 1, sum)
  ds$pred <- ifelse(ds$pred >= 2, 1,0)
  ds0<-data.frame(ds,ds0)
  
  
  
  
  
  colnames(ds1)<-c("A","B","C","pred","true")
  colnames(ds0)<-c("A","B","C","pred","true")
  dat<-rbind(ds0,ds1)
  
  
  
  ta<-table(dat$true,dat$pred)
  ta<-table(factor(dat$true,c("0","1"),levels=c("0","1")),factor(dat$pred,c("0","1"),levels=c("0","1")))
  
  truppv1<-rbeta(1,ta[4]+alpha_ppv,ta[3]+beta_ppv)
  truese1<-rbeta(1,ta[4]+alpha_se,ta[2]+beta_se)
  
  ppv <- cbind(ppv,truppv1)
  se <- cbind(se,truese1)
  
  tf<-2*truppv1*truese1/(truese1+truppv1)
  daa<-dat[,-5]
  daa<-daa[,-4]
  
  
  sj31.gibbs <- blca.gibbs(daa, 2, burn.in = 500, iter = 2000, 
                           alpha = c(1/2), beta =c(1/2),delta = 1/2,start.vals = c("prior"))
  
  gibbs <- as.mcmc(sj31.gibbs)
  dm <- data.frame(gibbs)
  
  
  d<-NULL
  c<-NULL
  b<-NULL
  nomal<-NULL
  no<-NULL
  a<-NULL
  for(l in 1:nrow(dm)){
    b<-(dm$ClassProb.2[l]
        *(dm$ItemProb.2.1[l]^dat$A*((1-dm$ItemProb.2.1[l])^(1-dat$A)))
        *(dm$ItemProb.2.2[l]^dat$B*((1-dm$ItemProb.2.2[l])^(1-dat$B)))
        *(dm$ItemProb.2.3[l]^dat$C*((1-dm$ItemProb.2.3[l])^(1-dat$C)))
        
    )
    a<-(dm$ClassProb.1[l]
        *(dm$ItemProb.1.1[l]^dat$A*((1-dm$ItemProb.1.1[l])^(1-dat$A)))
        *(dm$ItemProb.1.2[l]^dat$B*((1-dm$ItemProb.1.2[l])^(1-dat$B)))
        *(dm$ItemProb.1.3[l]^dat$C*((1-dm$ItemProb.1.3[l])^(1-dat$C)))
    )
    
    if (dm$ItemProb.2.1[l]+dm$ItemProb.2.2[l]+dm$ItemProb.2.3[l]>dm$ItemProb.1.1[l]+dm$ItemProb.1.2[l]+dm$ItemProb.1.3[l]){
      d<-b/(a+b)
    }
    else {
      d<-a/(a+b)
    }
    c <- cbind(c,d)
  }
  
  
  d<-NULL
  asa1<-NULL
  az<-NULL
  for(j in 1:2000){
    d<-rbinom(n, 1,c[,j])
    asa1<-cbind(asa1,d)
  }
  asa<-data.frame(asa1)
  
  
  prop<-c
  a<-NULL
  b<-NULL
  c<-NULL
  d<-NULL
  g<-NULL
  aa<-NULL
  ab<-NULL
  ac<-NULL
  ae<-NULL
  preppv<-NULL
  prese<-NULL
  appv<-NULL
  ase<-NULL
  af<-NULL
  ma<-NULL
  for(t in 1:2000){
    g<-table(factor(asa[,t],c("0","1"),levels=c("0","1")),factor(dat$pred,c("0","1"),levels=c("0","1")))
    a<-g[1]
    b<-g[2]
    c<-g[3]
    e<-g[4]
    preppv<-rbeta(1, e+alpha_ppv, c+beta_ppv)
    prese<-rbeta(1,e+alpha_se,b+beta_se)
    fs<-2*preppv*prese/(preppv+prese)
    af<-cbind(af,fs)
  }
  
  f1<-af
  bf<-sample(f1,1)
  (hdiMC <- hdi(as.vector(f1)))
  cov=c(cov,hdiMC[1]<tf&tf<hdiMC[2])
}



a<-summary(cov)
a1<-as.numeric(a[2])
a2<-as.numeric(a[3])
ds<-a2/(a1+a2)
b<-data.frame(ds,a1,a2)
colnames(b) <- c("cv","False","True")
b