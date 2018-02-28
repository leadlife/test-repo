setwd("D:/SUNY_SB-QF_PHD-Courses/Systemic risk credit risk/paperwork/ntest")
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
lambda <- seq(0,100,1);
NStart <- 75;
NRateY <- 36;
NFactEnd <- N;
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
 Xmatr <- cbind(hermite(tempx,0),hermite(tempx,1),hermite(tempx,2),hermite(tempx,3),hermite(tempx,4));
 Varyhat <- (1+Xmatr%*%CovBeta%*%t(Xmatr))*Covres[PVSelIndex[i]];
 Beta <- LambdaSel[[6]][[PVSelIndex[i]]][[1]];
 ##Tangentch <- TangentCheck1(Xmatr,tempx,Beta);
 ##yhat[i] <- Tangentch[[1]];
 yhat[i] <- Xmatr%*%Beta;
 ##Slope[i] <- Tangentch[[2]];
 ##Intercept[i] <- Tangentch[[3]];
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

Beta <- LambdaSel[[6]][[MinErrIndex]][[1]];
tempx <- -1;
Tangentch <- TangentCheck1(tempx,Beta);
InterceptN1 <- Tangentch[[3]];
SlopeN1 <- Tangentch[[2]];
tempx <- 0;
Tangentch0 <- TangentCheck1(tempx,Beta);
Intercept0 <- Tangentch0[[3]];
Slope0 <- Tangentch0[[2]];
MinErrIndex 
IR1max
yhat1	

## draw curve
plot(LambdaSel[[5]][[MinErrIndex]],LambdaSel[[4]][[MinErrIndex]],xlab="Best Factor",ylab="Citi(Normal period)");
datax <- LambdaSel[[5]][[MinErrIndex]];
datay <- LambdaSel[[6]][[MinErrIndex]][[6]]%*%LambdaSel[[6]][[MinErrIndex]][[1]];
##plot(datax,datay,type="n",xlab="X",ylab="Y");
drawlines(datax,datay);
ablineclip(a=InterceptN1,b=SlopeN1,x1=-1.3,x2=-0.7,col="red",lwd=2);
ablineclip(a=Intercept0,b=Slope0,x1=-0.3,x2=0.3,col="black",lwd=2);
legend("topleft",legend=c("Delta at -1","Delta at 0"),col=c("red","black"),lwd=2,lty=1, cex=0.7);

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







