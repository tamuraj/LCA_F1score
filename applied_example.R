#-----------------------------------------------------------------------------------------------------#
#   Applied Example: Alzheimer data in Moran et al (2004)               										          #
#-----------------------------------------------------------------------------------------------------#


library(HDInterval)
library(BayesLCA)
# Start measuring time
start_time <- Sys.time()
#----------------------------------------------------------------------------------------------------------------------#
#   Assessing the diagnostic performance of Alzheimer's disease by Hallucination by giving an estimate of the F1 score #				  
#----------------------------------------------------------------------------------------------------------------------#



# dat         dataframe       data set
# test	  vector		Hallucination' response



data(Alzheimer)

dat<-Alzheimer[,-1]

test<-Alzheimer[,1]



sj31.gibbs <- blca.gibbs(dat, 2, burn.in = 500, iter = 2000, 
                         alpha = c(1/2), beta =c(1/2),delta = 1/2,start.vals = c("prior"))
gibbs <- as.mcmc(sj31.gibbs)
dm <- data.frame(gibbs)

#plot(sj31.gibbs,3)


# Probability of belonging to Class 2
lcm1<- NULL
b<- NULL
for(i in 1:nrow(dm)){
  b<-(
    dm$ClassProb.2[i]
    *(dm$ItemProb.2.1[i]^dat$Activity*((1-dm$ItemProb.2.1[i])^(1-dat$Activity)))
    *(dm$ItemProb.2.2[i]^dat$Aggression*((1-dm$ItemProb.2.2[i])^(1-dat$Aggression)))
    *(dm$ItemProb.2.3[i]^dat$Agitation*((1-dm$ItemProb.2.3[i])^(1-dat$Agitation)))
    *(dm$ItemProb.2.4[i]^dat$Diurnal*((1-dm$ItemProb.2.4[i])^(1-dat$Diurnal)))
    *(dm$ItemProb.2.5[i]^dat$Affective*((1-dm$ItemProb.2.5[i])^(1-dat$Affective)))
  )/
    (
      (
        dm$ClassProb.1[i]
        *(dm$ItemProb.1.1[i]^dat$Activity*((1-dm$ItemProb.1.1[i])^(1-dat$Activity)))
        *(dm$ItemProb.1.2[i]^dat$Aggression*((1-dm$ItemProb.1.2[i])^(1-dat$Aggression)))
        *(dm$ItemProb.1.3[i]^dat$Agitation*((1-dm$ItemProb.1.3[i])^(1-dat$Agitation)))
        *(dm$ItemProb.1.4[i]^dat$Diurnal*((1-dm$ItemProb.1.4[i])^(1-dat$Diurnal)))
        *(dm$ItemProb.1.5[i]^dat$Affective*((1-dm$ItemProb.1.5[i])^(1-dat$Affective)))
      )
      +(
        dm$ClassProb.2[i]
        *(dm$ItemProb.2.1[i]^dat$Activity*((1-dm$ItemProb.2.1[i])^(1-dat$Activity)))
        *(dm$ItemProb.2.2[i]^dat$Aggression*((1-dm$ItemProb.2.2[i])^(1-dat$Aggression)))
        *(dm$ItemProb.2.3[i]^dat$Agitation*((1-dm$ItemProb.2.3[i])^(1-dat$Agitation)))
        *(dm$ItemProb.2.4[i]^dat$Diurnal*((1-dm$ItemProb.2.4[i])^(1-dat$Diurnal)))
        *(dm$ItemProb.2.5[i]^dat$Affective*((1-dm$ItemProb.2.5[i])^(1-dat$Affective)))
      )
    )
  
  lcm1 <- cbind(lcm1,b)
  
}



lcmdat<-data.frame(lcm1)

d<-NULL
az<-NULL
ds<-NULL

for(i in 1:2000){
  dst<-lcm1[,i]
  d<-rbinom(240, 1,dst)
  ds<-cbind(ds,d)
}

lcm2<-matrix(ds,nrow = 240,ncol = 2000)
lcmdat<-data.frame(lcm2)






ppv <- se <- sp <- npv <- ac <- numeric(2000)

for (i in 1:2000) {
  dst <- lcm1[, i]
  pred <- rbinom(240, 1, dst)
  
  g <- table(factor(test, levels = 0:1), factor(pred, levels = 0:1))
  TP <- g["1", "1"]
  FN <- g["0", "1"]
  FP <- g["1", "0"]
  TN <- g["0", "0"]
  
  # Jeffreys priorによるベータ分布のサンプリング
  ppv[i] <- rbeta(1, TP + 0.5, FP + 0.5)
  se[i]  <- rbeta(1, TP + 0.5, FN + 0.5)
  sp[i]  <- rbeta(1, TN + 0.5, FP + 0.5)
  npv[i] <- rbeta(1, TN + 0.5, FN + 0.5)
  ac[i]  <- (TP + TN) / (TP + TN + FP + FN)
}

f1<-2*ppv*se/(ppv+se)


(hdiMC <- hdi(as.vector(f1)))


plot(density(f1), xlim=c(0,1), col = "red", lwd = 2,xlab="",main="",ylab="")
abline(v=hdiMC,col="blue")  ##95% Highest density interval


# End measuring time
end_time <- Sys.time()

# Calculate and print the elapsed time
execution_time <- end_time - start_time
print(paste("Execution time:", execution_time))
