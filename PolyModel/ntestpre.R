setwd("D:/1-SUNY_SB-QF_PHD-Courses/Systemic risk credit risk/paperwork/ntest")
library(zoo)
library(EQL)
library(MASS)
library(base)
##library(broom)
library(matrixcalc)
library(xts)
library(pracma)
library(polynom)
library(orthopolynom)
library(numDeriv)
library(stats)
library(plotrix)
library()

source("D:/1-SUNY_SB-QF_PHD-Courses/Systemic risk credit risk/paperwork/ntest/FacDtrnVaR.R");
source("D:/1-SUNY_SB-QF_PHD-Courses/Systemic risk credit risk/paperwork/ntest/OneFacRidgReg.R");
source("D:/1-SUNY_SB-QF_PHD-Courses/Systemic risk credit risk/paperwork/ntest/drawlines.R");
source("D:/1-SUNY_SB-QF_PHD-Courses/Systemic risk credit risk/paperwork/ntest/DataPro.R");
source("D:/1-SUNY_SB-QF_PHD-Courses/Systemic risk credit risk/paperwork/ntest/MulFacRidReg.R");
source("D:/1-SUNY_SB-QF_PHD-Courses/Systemic risk credit risk/paperwork/ntest/Merge2RidgReg.R");
source("D:/1-SUNY_SB-QF_PHD-Courses/Systemic risk credit risk/paperwork/ntest/TangentCheck1.R");
source("D:/1-SUNY_SB-QF_PHD-Courses/Systemic risk credit risk/paperwork/ntest/TangentCheck2.R");
source("D:/1-SUNY_SB-QF_PHD-Courses/Systemic risk credit risk/paperwork/ntest/zmatrixpro.R");

# Data Loading
FIndex <- read.csv('Factor Index518.csv', header = TRUE, stringsAsFactors = FALSE)
N <- length(FIndex$Symbol)
factor_data <- list()
respond_data <- list()
Len <- rep(0,N);
FIndex_file <- vector(mode = "character", length = N)
#set.seed(2);
for (i in 1:N)
{
  FIndex_file[i] <- paste(FIndex$Symbol[i], '.csv', sep = "")
  factor_data[[i]] <- read.csv(FIndex_file[i], header = TRUE, stringsAsFactors = FALSE, skip = 1)
  if(i>=151)
  {
   factor_data[[i]]$Date <- as.Date(factor_data[[i]]$Date, format = "%m/%d/%Y")
  }
  factor_data[[i]]$Date <- as.Date(factor_data[[i]]$Date, format = "%m/%d/%y")
  names(factor_data[[i]])[2] <- FIndex$Symbol[i]
  Len[i] <- length(factor_data[[i]]$Date)
  
  if (FIndex$Frequency[i] == "M") 
  {
    factor_data[[i]] <- data.frame(factor_data[[i]], row.names = 1)
    factor_data[[i]] <- as.xts(factor_data[[i]])
    index(factor_data[[i]]) <- as.yearmon(index(factor_data[[i]]))
  } 
  else if (FIndex$Frequency[i] == "Q") 
  {
    factor_data[[i]] <- data.frame(factor_data[[i]], row.names = 1)
    factor_data[[i]] <- as.xts(factor_data[[i]])
    index(factor_data[[i]]) <- as.yearmon(index(factor_data[[i]]))
  } 
  else 
  {
    factor_data[[i]] <- NA
  }  
}

##one-factor hermite poly regression
hermitede <- function(x){cbind(0,1,2*x,-3+3*x^2,-12*x + 4*x^3)};
result <- list();
lambda <- seq(0,100,1);
NStart <- 100;
NRateY <- 36;
NFactEnd <- 179;
NFactBegin <- 1;
##Dependent data process
IndexDep <- 153; 
set.seed(1);
LambdaSel <- 0;
##Factor data process
Penamatrix <- matrix(c(0,0,0,0,0, 0,0,0,0,0, 0,0,4,0,0, 0,0,0,9,0, 0,0,0,0,16),5,5);
LambdaSel <- MulFacRidReg(lambda,factor_data,IndexDep,NFactBegin,NFactEnd,NStart,NRateY);
ErrDiff <- LambdaSel[[2]];
##ErrRatio <- LambdaSel[[3]];
LSErr <- LambdaSel[[1]];
RSquared <- LambdaSel[[8]];
VPvalue <- LambdaSel[[7]];
Covres <- LambdaSel[[9]];
##PVSelIndex <- which(VPvalue <= sqrt(min(VPvalue[-c(IndexDep)])));
PVSelIndex <- which(VPvalue <= max(head(sort(VPvalue),20)));
##IR
NPVIndex <- length(PVSelIndex);
StValue <- -1;
IR <- rep(0,NPVIndex);
yhat <- rep(0,NPVIndex);
for(i in 1:NPVIndex)
{
 tempy <- LambdaSel[[4]][[PVSelIndex[i]]];
 tempx <- StValue;
 zmatrix <- LambdaSel[[6]][[PVSelIndex[i]]][[6]];
 CovBetatemp <- matrix.inverse(t(zmatrix)%*%zmatrix+LambdaSel[[6]][[PVSelIndex[i]]][[2]]*Penamatrix)%*%t(zmatrix);
 CovBeta <- CovBetatemp%*%t(CovBetatemp);
 Xmatr <- zmatrixpro(tempx);
 Varyhat <- (1+Xmatr%*%CovBeta%*%t(Xmatr))*Covres[PVSelIndex[i]];
 Beta <- LambdaSel[[6]][[PVSelIndex[i]]][[1]];
 yhat[i] <- Xmatr%*%Beta;
 IR[i] <- sqrt(yhat[i]^2/Varyhat);
}

IR[which(PVSelIndex==IndexDep)] <- 0;
IRMaxIndex <- which.max(IR);
IR1max <- max(IR);
MinErrIndex <- PVSelIndex[IRMaxIndex];
yhat1 <- yhat[IRMaxIndex];
IRSelIndex <- which(IR >= min(tail(sort(IR),10)));
IRorder <- rev(tail(sort(IR),10));
NIRorder <- length(IRorder)
IRorderIndex <- rep(0,NIRorder);
NFacdata <- list();
for( i in 1:NIRorder)
{
 IRorderIndex[i] <- which(IR == IRorder[i]);
 temp <- PVSelIndex[IRorderIndex[i]];
 NFacdata[[i]] <- LambdaSel[[5]][[temp]];
}
MinErrIndex
PVSelIndex[IRorderIndex]
NDepdata <- LambdaSel[[4]][[MinErrIndex]];
Penamatrix2 <- matrix(0, nrow = 10, ncol = 10);
for (i in 1:2)
{
 for (k in 1:5)
 {
  Penamatrix2[5*(i-1)+k,5*(i-1)+k] <- (k-1)^2;
 }
}
Penamatrix2[2,2] <- 0;
Penamatrix2 <- Penamatrix2[-6,-6];
Penamatrix2[6,6] <- 0;
IR2 <- matrix(0, nrow = NIRorder-1, ncol = NIRorder);
yhat2 <- matrix(0, nrow = NIRorder-1, ncol = NIRorder);
for( i in 1:(NIRorder-1))
{
 NFacdata1 <- NFacdata[[i]];
 for (j in (i+1):NIRorder)
 {
  NFacdata2 <- NFacdata[[j]];
  Merge2Result <- Merge2RidgReg(NDepdata,NFacdata1,NFacdata2,lambda);
  tempy <- NDepdata;
  tempx1 <- StValue;
  tempx2 <- StValue;
  zmatrix <- Merge2Result[[6]];
  CovBetatemp <- matrix.inverse(t(zmatrix)%*%zmatrix+Merge2Result[[2]]*Penamatrix2)%*%t(zmatrix);
  CovBeta <- CovBetatemp%*%t(CovBetatemp);
  Xmatr1 <- zmatrixpro(tempx1);
  Xmatr2 <- zmatrixpro(tempx2);
  Xmatr <- matrix(cbind(Xmatr1,Xmatr2)[1,-6]);
  Varyhat2 <- (1+t(Xmatr)%*%CovBeta%*%Xmatr)*Merge2Result[[8]];
  yhat2[i,j] <- t(Xmatr)%*%Merge2Result[[1]];
  IR2[i,j] <- sqrt(yhat2[i,j]^2/Varyhat2);
 }
}
MinErrIndex 
IR1max
yhat1	
IR2maxloc <- which(IR2==max(IR2),arr.ind = TRUE);
IR2maxloc
FinafacIndex <- PVSelIndex[IRorderIndex[IR2maxloc]];
FinafacIndex
max(IR2)
yhat2[IR2maxloc]
## draw curve
plot(LambdaSel[[5]][[MinErrIndex]],LambdaSel[[4]][[MinErrIndex]]);
datax <- LambdaSel[[5]][[MinErrIndex]];
datay <- LambdaSel[[6]][[MinErrIndex]][[6]]%*%LambdaSel[[6]][[MinErrIndex]][[1]];
drawlines(datax,datay);
x11();
plot(LambdaSel[[5]][[FinafacIndex[1]]],LambdaSel[[4]][[FinafacIndex[1]]]);
datax <- LambdaSel[[5]][[FinafacIndex[1]]];
datay <- LambdaSel[[6]][[FinafacIndex[1]]][[6]]%*%LambdaSel[[6]][[FinafacIndex[1]]][[1]];
drawlines(datax,datay);
x11();
plot(LambdaSel[[5]][[FinafacIndex[2]]],LambdaSel[[4]][[FinafacIndex[2]]]);
datax <- LambdaSel[[5]][[FinafacIndex[2]]];
datay <- LambdaSel[[6]][[FinafacIndex[2]]][[6]]%*%LambdaSel[[6]][[FinafacIndex[2]]][[1]];
drawlines(datax,datay);
##VaR calculation: 5 points(99%,84%,50%,16%,1%)
source("//mysbfiles.campus.stonybrook.edu/~/ntest/FacDtrnVaR.R");
Ximin <- 1.5; #1.5*sigma
nVaR <- list();

for (i in 1:N)
{
 nVaR[[i]] <- FacDtrnVaR(factor_data[[i]],Ximin,i,Len[i]-1);
 ##colnames(nVaR[[i]]) <- FIndex$Symbol[i];
} 
nVaR
##save(result,Errvector,MinErrIndex,Bestlambda,file = "result.RData");
##load("result.RData")







