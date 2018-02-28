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

source("D:/SUNY_SB-QF_PHD-Courses/Systemic risk credit risk/paperwork/ntest/FacDtrnVaR.R");
source("D:/SUNY_SB-QF_PHD-Courses/Systemic risk credit risk/paperwork/ntest/OneFacRidgReg.R");
source("D:/SUNY_SB-QF_PHD-Courses/Systemic risk credit risk/paperwork/ntest/drawlines.R");
source("D:/SUNY_SB-QF_PHD-Courses/Systemic risk credit risk/paperwork/ntest/DataPro.R");
source("D:/SUNY_SB-QF_PHD-Courses/Systemic risk credit risk/paperwork/ntest/MulFacRidReg.R");
source("D:/SUNY_SB-QF_PHD-Courses/Systemic risk credit risk/paperwork/ntest/Merge2RidgReg.R");
source("D:/SUNY_SB-QF_PHD-Courses/Systemic risk credit risk/paperwork/ntest/TangentCheck1.R");
source("D:/SUNY_SB-QF_PHD-Courses/Systemic risk credit risk/paperwork/ntest/TangentCheck2.R");
source("D:/SUNY_SB-QF_PHD-Courses/Systemic risk credit risk/paperwork/ntest/zmatrixpro.R");

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
hermitede <- function(x){c(1,x,2*x,-2,-8*x + 6)};
result <- list();
NP <- 121;
lambda <- seq(0,100,1);
Tan <- matrix(0,NP,3);
Vyhat <- matrix(0,NP,1);
NRateY <- 30;
NFactEnd <- 179;
NFactBegin <- 1;
##Dependent data process 

VTanPN1 <- matrix(0,NP,179-152+1);
VTanP0 <- matrix(0,NP,179-152+1);
VTanPP1 <- matrix(0,NP,179-152+1);
VVyhat <- matrix(0,NP,179-152+1);
for(m in 1:28)
{
  IndexDep <- 151+m;
  print(IndexDep);
  for(j in 1:NP)
  {
    NStart <- 122-j;
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
    VPvalue[PVSelIndex];
    ##IR
    NPVIndex <- length(PVSelIndex);
    StValue <- -1;
    IR <- rep(0,NPVIndex);
    yhat <- rep(0,NPVIndex);
    slope <- rep(0,NPVIndex);
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
    IR1max
    Vyhat[j] <- yhat1;	
    BetaM <- LambdaSel[[6]][[MinErrIndex]][[1]]
    deltax <- 0.0005;
    TanPN1 <- (zmatrixpro(-1+deltax)%*%BetaM-zmatrixpro(-1-deltax)%*%BetaM)/2/deltax;
    TanP0 <- (zmatrixpro(0+deltax)%*%BetaM-zmatrixpro(0-deltax)%*%BetaM)/2/deltax;
    TanPP1 <- (zmatrixpro(1+deltax)%*%BetaM-zmatrixpro(1-deltax)%*%BetaM)/2/deltax;
    Tan[j,] <- c(TanPN1,TanP0,TanPP1);
  }
  VTanPN1[,m] <- Tan[,1];
  VTanP0[,m] <- Tan[,2];
  VTanPP1[,m] <- Tan[,3];
  VVyhat[,m] <- Vyhat;
}
VTanPN1C <- VTanPN1;
VTanPN1C[which(VTanPN1C<0)] <- abs(VTanPP1[which(VTanPN1C<0)]);
VTanP0C <- abs(VTanP0);
Date <- rev(rev(index(factor_data[[1]]))[1:119]);
plot(cbind(Date,VTanPN1C[,2]),type="l",col="red",lwd=2,xlab="Month-Year",ylab="AD",xaxt = "n");
axis(1,at = Date,labels = strftime(Date, format="%m/%Y"),cex.axis=0.5,las=2);
lines(Date,VTanPN1C[,5],pch=".",col="blue", lwd=2,type="l", lty=1);
legend("topright", legend=c("Citi Delta at -1","Morgan Stanley Delta at -1"),col=c("red","blue"),lwd=2,lty=1, cex=0.6);
x11();
plot(Date,VTanP0C[,2],type="l",col="red",lwd=2,xlab="Month-Year",ylab="AD",xaxt = "n");
axis(1,at = Date,labels = strftime(Date, format="%m/%Y"),cex.axis=0.5,las=2);
lines(Date,VTanP0C[,5],pch=".",col="blue", lwd=2,type="l", lty=1);
legend("topright", legend=c("Citi Delta at 0","Morgan Stanley Delta at 0"),col=c("red","blue"),lwd=2,lty=1, cex=0.6);
x11();
ADPN1 <- rowSums(VTanPN1C)/20;
plot(Date,ADPN1,type="l",col="red",lwd=2,xlab="Month-Year",ylab="AD at -1*sigma",xaxt = "n");
axis(1,at = Date,labels = strftime(Date, format="%m/%Y"),cex.axis=0.5,las=2);
legend("topright", legend="AD at -1",col="red", lwd=2,lty=1, cex=0.8);
x11();
ADP0 <- rowSums(VTanP0C)/20;
plot(Date,ADP0,type="l",col="blue",lwd=2,xlab="Month-Year",ylab="AD at 0",xaxt = "n");
axis(1,at = Date,labels = strftime(Date, format="%m/%Y"),cex.axis=0.5,las=2);
legend("topright", legend="AD at 0",col="blue", lwd=2,lty=1, cex=0.8);
x11();
ASV <- -1.5*rowSums(abs(VVyhat))/20;
plot(Date,ASV,type="l",col="blue",lwd=2,xlab="Month-Year",ylab="ASV",xaxt = "n");
axis(1,at = Date,labels = strftime(Date, format="%m/%Y"),cex.axis=0.5,las=2);
legend("top", legend="ASV",col="blue", lwd=2,lty=1, cex=0.8);
save(VTanPN1,VTanP0,VTanPP1,VVyhat,file = "VTan.RData");
load("VTan.RData")

CSV <- -abs(VVyhat[,5]);
plot(Date,CSV,type="l",col="blue",lwd=2,xlab="Month-Year",ylab="Morgan Stanley Stress VaR",xaxt = "n");
axis(1,at = Date,labels = strftime(Date, format="%m/%Y"),cex.axis=0.5,las=2);
legend("top", legend="Morgan Stanley Stress VaR",col="blue", lwd=2,lty=1, cex=0.6);
## draw curve
plot(LambdaSel[[5]][[MinErrIndex]],LambdaSel[[4]][[MinErrIndex]]);
datax <- LambdaSel[[5]][[MinErrIndex]];
datay <- LambdaSel[[6]][[MinErrIndex]][[6]]%*%LambdaSel[[6]][[MinErrIndex]][[1]];
drawlines(datax,datay);






