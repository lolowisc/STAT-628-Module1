#read data
dat <- read.csv("BodyFat.csv")
dat <- dat[,-1]

#data cleaning
#check siri's equation
plot(dat$BODYFAT, (495/dat$DENSITY-450), cex=0.1, xlim = c(-5,50), 
     xlab="bodyfat (%)", ylab="495 / density - 450",main="siri's equation")
text(dat$BODYFAT[c(48,76,96,182)],(495/dat$DENSITY[c(48,76,96,182)]-450),
     labels=c(48,76,96,182),cex=0.8,pos=4)
#check BMI equation
plot(dat$ADIPOSITY,(0.45359237*dat$WEIGHT)/(0.0254*dat$HEIGHT)^2,cex=0.1,
     xlab="adiposity",ylab="weight / height^2",main="bmi equation")
text(dat$ADIPOSITY[c(42,163,221)],(0.45359237*dat$WEIGHT[c(42,163,221)])/(0.0254*dat$HEIGHT[c(42,163,221)])^2,
     labels=c(42,163,221),pos=2,cex=0.8)
dat$HEIGHT[42] <- sqrt((0.45359237*dat$WEIGHT[42])/dat$ADIPOSITY[42])/0.0254 #fix the point
#check cook's distance
bodyfat <- dat[,-c(2)]
summary(model <- lm(BODYFAT ~ ., data=bodyfat[-c(48,76,96,182),]))
layout(matrix(1:1, ncol=1))
plot(model, which=4)
abline( h = 4/(248-15),lty=2 )
#after deleting outliers
summary(model <- lm(BODYFAT ~ ., data=bodyfat[-c(48,76,96,182,39,86,221),]))
layout(matrix(1:1, ncol=1))
plot(model, which=4)
abline( h = 4/(245-15),lty=2 )

#data set after cleaning
bodyfat <- bodyfat[-c(48,76,96,182,39,86,221),]; 
rownames(bodyfat) <- NULL

result <- data.frame(matrix(NA,nrow=10,ncol=5))
rownames(result) <- c("original","cp","r2","back.AIC","back.BIC","forw.AIC","forw.BIC","rm","lasso","bess")
colnames(result) <- c("cv.mse","mse","ad.r2","aic","vif")

#####
#check mse through cross validation
library(DAAG)
library(car)
#original model with full variables
model <- lm(BODYFAT ~ ., data=bodyfat)
model_cv <- cv.lm(BODYFAT ~ AGE + WEIGHT + HEIGHT + ADIPOSITY + NECK + CHEST + 
                    ABDOMEN + HIP + THIGH + KNEE + ANKLE + BICEPS + FOREARM + 
                    WRIST , data=bodyfat, m=10, printit = F) 
result$cv.mse[1] <- mean((bodyfat$BODYFAT-model_cv$cvpred)^2)
result$mse[1] <- mean((model$residuals)^2)
result$ad.r2[1] <- summary(model)$adj.r.squared
result$aic[1] <- extractAIC(model)[2]
result$vif[1] <- paste(paste(names(sort(vif(model),decreasing = T)),
                              round(sort(vif(model),decreasing = T),2),sep=":"),collapse = " ")

#####
#Mallow's Cp
X <- model.matrix(model)[,-1]
Y <- bodyfat[,1]
library(leaps) # for leaps()
library(faraway) # for Cpplot()
g <- leaps(X, Y, nbest=1)
layout(matrix(1:1, ncol=1))
Cpplot(g)
cp.choice <- c(1,2,7,9,13,14)+1
summary(model.cp <- lm(BODYFAT ~ ., data=bodyfat[,c(1,cp.choice)]))
vif(model.cp)
form <- as.formula(paste("BODYFAT~",
                         paste(colnames(bodyfat)[cp.choice],collapse = "+"),sep=""))
cp_cv <- cv.lm(form, data=bodyfat, m=10, printit = F)
result$cv.mse[2] <- mean((bodyfat$BODYFAT-cp_cv$cvpred)^2)
result$mse[2] <- mean((model.cp$residuals)^2)
result$ad.r2[2] <- summary(model.cp)$adj.r.squared
result$aic[2] <- extractAIC(model.cp)[2]
result$vif[2] <- paste(paste(names(sort(vif(model.cp),decreasing = T)),
                             round(sort(vif(model.cp),decreasing = T),2),sep=":"),collapse = " ")


#####
#r^2
g <- leaps(X,Y, nbest=1, method="adjr2")
plot(g$adjr2)
(g$which)[which(g$adjr2 == max(g$adjr2)),]
r2.choice <- which((g$which)[which(g$adjr2 == max(g$adjr2)),])+1
form <- as.formula(paste("BODYFAT~",
                         paste(colnames(bodyfat)[r2.choice],collapse = "+"),sep=""))
summary(model.r2 <- lm(BODYFAT ~ ., data=bodyfat[,c(1,r2.choice)]))
r2_cv <- cv.lm(form, data=bodyfat, m=10, printit = F)
result$cv.mse[3] <- mean((bodyfat$BODYFAT-r2_cv$cvpred)^2)
result$mse[3] <- mean((model.r2$residuals)^2)
result$ad.r2[3] <- summary(model.r2)$adj.r.squared
result$aic[3] <- extractAIC(model.r2)[2]
result$vif[3] <- paste(paste(names(sort(vif(model.r2),decreasing = T)),
                             round(sort(vif(model.r2),decreasing = T),2),sep=":"),collapse = " ")

#####
#backward selection
summary(model.AIC <- step(model,direction = "both" ,k=2,trace = 0))
summary(model.BIC <- step(model,direction = "both", k=log(245),trace = 0))
AIC.choice <- which(colnames(bodyfat)%in%names(model.AIC$coefficients)[-1])
form <- as.formula(paste("BODYFAT~",
                         paste(colnames(bodyfat)[AIC.choice],collapse = "+"),sep=""))
AIC_cv <- cv.lm(form, data=bodyfat, m=10, printit = F)
result$cv.mse[4] <- mean((bodyfat$BODYFAT-AIC_cv$cvpred)^2)
result$mse[4] <- mean((model.AIC$residuals)^2)
result$ad.r2[4] <- summary(model.AIC)$adj.r.squared
result$aic[4] <- extractAIC(model.AIC)[2]
result$vif[4] <- paste(paste(names(sort(vif(model.AIC),decreasing = T)),
                             round(sort(vif(model.AIC),decreasing = T),2),sep=":"),collapse = " ")

BIC.choice <- which(colnames(bodyfat)%in%names(model.BIC$coefficients)[-1])
form <- as.formula(paste("BODYFAT~",
                         paste(colnames(bodyfat)[BIC.choice],collapse = "+"),sep=""))
BIC_cv <- cv.lm(form, data=bodyfat, m=10, printit = F)
result$cv.mse[5] <- mean((bodyfat$BODYFAT-BIC_cv$cvpred)^2)
result$mse[5] <- mean((model.BIC$residuals)^2)
result$ad.r2[5] <- summary(model.BIC)$adj.r.squared
result$aic[5] <- extractAIC(model.BIC)[2]
result$vif[5] <- paste(paste(names(sort(vif(model.BIC),decreasing = T)),
                             round(sort(vif(model.BIC),decreasing = T),2),sep=":"),collapse = " ")

#forward selection
base <- lm(BODYFAT~1,data=bodyfat) 
summary(AIC.base<-step(base,direction="both", 
               scope=list(lower=~1,upper=model),trace=0))
summary(BIC.base<-step(base,direction="both", k=log(245),
                       scope=list(lower=~1,upper=model),trace=0))
AICbs.choice <- which(colnames(bodyfat)%in%names(AIC.base$coefficients)[-1])
form <- as.formula(paste("BODYFAT~",
                         paste(colnames(bodyfat)[AICbs.choice],collapse = "+"),sep=""))
AICbs_cv <- cv.lm(form, data=bodyfat, m=10, printit = F)
result$cv.mse[6] <- mean((bodyfat$BODYFAT-AICbs_cv$cvpred)^2)
result$mse[6] <- mean((AIC.base$residuals)^2)
result$ad.r2[6] <- summary(AIC.base)$adj.r.squared
result$aic[6] <- extractAIC(AIC.base)[2]
result$vif[6] <- paste(paste(names(sort(vif(AIC.base),decreasing = T)),
                             round(sort(vif(AIC.base),decreasing = T),2),sep=":"),collapse = " ")

BICbs.choice <- which(colnames(bodyfat)%in%names(BIC.base$coefficients)[-1])
form <- as.formula(paste("BODYFAT~",
                         paste(colnames(bodyfat)[BICbs.choice],collapse = "+"),sep=""))
BICbs_cv <- cv.lm(form, data=bodyfat, m=10, printit = F)
result$cv.mse[7] <- mean((bodyfat$BODYFAT-BICbs_cv$cvpred)^2)
result$mse[7] <- mean((BIC.base$residuals)^2)
result$ad.r2[7] <- summary(BIC.base)$adj.r.squared
result$aic[7] <- extractAIC(BIC.base)[2]
summary(lm(form,data=bodyfat))
result$vif[7] <- paste(paste(names(sort(vif(BIC.base),decreasing = T)),
                             round(sort(vif(BIC.base),decreasing = T),2),sep=":"),collapse = " ")

#####
#random forest
#install.packages("party")
library(party)
set.seed(2019)
crf<-cforest(BODYFAT~.,control = cforest_unbiased(mtry = 2, ntree = 50), data=bodyfat)  
varimpt<-data.frame(varimp(crf))  
plot(sort(varimpt[,1],decreasing = T))
variables <- rownames(varimpt)[order(varimpt,decreasing = T)[1:6]]
form <- as.formula(paste("BODYFAT~",
                         paste(variables,collapse = "+"),sep=""))
summary(model.rm <- lm(form,data=bodyfat))
rm_cv <- cv.lm(form,data=bodyfat, m=10, printit = F)
result$cv.mse[8] <- mean((bodyfat$BODYFAT-rm_cv$cvpred)^2)
result$mse[8] <- mean((model.rm$residuals)^2)
result$ad.r2[8] <- summary(model.rm)$adj.r.squared
result$aic[8] <- extractAIC(model.rm)[2]
result$vif[8] <- paste(paste(names(sort(vif(model.rm),decreasing = T)),
                             round(sort(vif(model.rm),decreasing = T),2),sep=":"),collapse = " ")

#####
#LASSO
set.seed(2019)
library(glmnet)
cv.lasso <- cv.glmnet(x=data.matrix(bodyfat[,-1]), y=as.vector(bodyfat[,1]),nfold=10, alpha=1)
# Results
plot(cv.lasso)
plot(cv.lasso$glmnet.fit, xvar="lambda", label=TRUE)
cv.lasso$lambda.min
cv.lasso$lambda.1se
m_lasso=as.matrix(coef(cv.lasso, s=cv.lasso$lambda.1se))
m_lasso[!m_lasso[,1]==0,]
lasso_choice<-as.numeric(which(!m_lasso[,1]==0))[-1]
form <- as.formula(paste("BODYFAT~",
                         paste(colnames(bodyfat)[lasso_choice],collapse = "+"),sep=""))
summary(m_lasso <- lm(form,data=bodyfat))
par(mfrow=c(2,2));plot(m_lasso)
lasso_cv <- cv.lm(form,data=bodyfat, m=10, printit = F)
result$cv.mse[9] <- mean((bodyfat$BODYFAT-lasso_cv$cvpred)^2)
result$mse[9] <- mean((m_lasso$residuals)^2)
result$ad.r2[9] <- summary(m_lasso)$adj.r.squared
result$aic[9] <- extractAIC(m_lasso)[2]
result$vif[9] <- paste(paste(names(sort(vif(m_lasso),decreasing = T)),
                             round(sort(vif(m_lasso),decreasing = T),2),sep=":"),collapse = " ")

#BESS
#install.packages("BeSS")
par(mfrow=c(1,1))
library(BeSS)
bic.bess<-c()
for(i in 1:13){
  bess<-bess.one(x=data.matrix(bodyfat[,-1]), y=as.vector(bodyfat[,1]),s=i)
  bic.bess[i]<-bess$BIC
}
which(bic.bess==min(bic.bess))
bess<-bess.one(x=data.matrix(bodyfat[,-1]), y=as.vector(bodyfat[,1]),s=5)
summary(bess)
summary(m_bess <- lm(data=bodyfat,BODYFAT ~ ADIPOSITY + CHEST + ABDOMEN + WEIGHT + WRIST ))
bess_cv <- cv.lm(BODYFAT ~ ADIPOSITY + CHEST + ABDOMEN + WEIGHT + WRIST ,data=bodyfat, m=10, printit = F)
result$cv.mse[10] <- mean((bodyfat$BODYFAT-bess_cv$cvpred)^2)
result$mse[10] <- mean((m_bess$residuals)^2)
result$ad.r2[10] <- summary(m_bess)$adj.r.squared
result$aic[10] <- extractAIC(m_bess)[2]
result$vif[10] <- paste(paste(names(sort(vif(m_bess),decreasing = T)),
                             round(sort(vif(m_bess),decreasing = T),2),sep=":"),collapse = " ")

write.csv(result,file="comparison.csv")

#visualize the comparison between models
com <- result
com <- read.csv("comparison.csv")[,-c(3,5)]
colnames(com)[1] <- "method"
par(mar=c(8,4,4,5) + 0.1)
## Plot first set of data and draw its axis
plot(1:10, com$ad.r2, pch=16, axes=FALSE, ylim=c(0.7,0.8), xlab="", ylab="", type="b",col="black", main="Model Comparison")
axis(2, ylim=c(0.7,0.8),col="black",las=1)  ## las=1 makes horizontal labels
mtext("Adjusted R sqaure",side=2,line=3)
box()
par(new=TRUE)
## Plot the second plot and put axis scale on right
options(repr.plot.width = 5, repr.plot.height = 4)
plot(1:10, com$cv.mse, pch=15,  xlab="", ylab="", ylim=c(15,17), 
     axes=FALSE, type="b", col="red")
## a little farther out (line=4) to make room for labels
mtext("CV.MSE",side=4,col="red",line=3) 
axis(4, ylim=c(15,17), col="red",col.axis="red",las=1)
axis(1,1:10,labels = F)
text(1:10, par("usr")[3] - 0.25, srt = 45, adj = 1,
     labels = as.character(com$method), xpd = TRUE)
mtext("Models",side=1,col="black",line=6)

#reveal the relationship between wrist and bodyfat given abdomen to be almost the same
bodyfat$abdomen_level <- 1
bodyfat$abdomen_level[bodyfat$ABDOMEN>85.2] <- 2
bodyfat$abdomen_level[bodyfat$ABDOMEN>91] <- 3
bodyfat$abdomen_level[bodyfat$ABDOMEN>99.7] <- 4
ab1 <- bodyfat[bodyfat$abdomen_level==1,]
p1 <- ggplot(ab1, aes(WRIST, BODYFAT))+geom_smooth(method="lm")+
  labs(title ="abdomen in [70.4,85.2] cm", x = "wrist (cm)", y = "bodyfat (%)") + theme_classic()
ab2 <- bodyfat[bodyfat$abdomen_level==2,]
p2 <- ggplot(ab2, aes(WRIST, BODYFAT))+geom_smooth(method="lm")+
  labs(title ="abdomen in (85.2,91] cm", x = "wrist (cm)", y = "bodyfat (%)") + theme_classic()
ab3 <- bodyfat[bodyfat$abdomen_level==3,]
p3 <- ggplot(ab3, aes(WRIST, BODYFAT))+geom_smooth(method="lm")+
  labs(title ="abdomen in (91,99.7] cm", x = "wrist (cm)", y = "bodyfat (%)") + theme_classic()
ab4 <- bodyfat[bodyfat$abdomen_level==4,]
p4 <- ggplot(ab4, aes(WRIST, BODYFAT))+geom_smooth(method="lm")+
  labs(title ="abdomen in (99.7,126.2] cm", x = "wrist (cm)", y = "bodyfat (%)") + theme_classic()
library(gridExtra)
grid.arrange(p1, p2,p3,p4, nrow = 1)
