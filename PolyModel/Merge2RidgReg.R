Merge2RidgReg <- function(NDepdata,NFacdata1,NFacdata2,lambda)
{
 y <- NDepdata;
 Nsample <- length(y);
 LSerr <- rep(0,length(lambda));
 Penaltyerr <- rep(0,length(lambda));
 Differr <- rep(0,length(lambda));
 Beta <- list();
 tempx1 <- NFacdata1;
 tempx2 <- NFacdata2;
 zmatrix1 <- zmatrixpro(tempx1);
 tempy <- NDepdata;
 zmatrix2 <- zmatrixpro(tempx2);
 zmatrix <- cbind(zmatrix1,zmatrix2)[,-6];
 for (i in 1:length(lambda))
 {
   S <- matrix.inverse(t(zmatrix)%*%zmatrix+lambda[i]*Penamatrix2)%*%t(zmatrix);
   Beta[[i]] <- S%*%tempy;
   Fittedvalue <- zmatrix%*%Beta[[i]];
   LSerr[i] <- sum((Fittedvalue - tempy)^2);
   Penaltyerr[i] <- lambda[i]*diag(Penamatrix2)^2%*%Beta[[i]]^2; 
   
 }
 ratio <- abs(Penaltyerr/LSerr-1);
 temp <- which(ratio <= max(head(sort(ratio),10)));
 lambdaIndex <- temp[which.min(LSerr[temp])];
 Bestlambda <- lambda[lambdaIndex];
 SSres <- sum((zmatrix%*%Beta[[lambdaIndex]]- tempy)^2);
 covres <- SSres/(length(y)-5);
 SStot <- sum((y-mean(y))^2);
 RSquare <- 1-SSres/SStot;
 m <- 0;
 S <- matrix.inverse(t(zmatrix)%*%zmatrix+Bestlambda*Penamatrix2)%*%t(zmatrix);

for (j in 1:1000)
 {
  PermuY <- sample(tempy,NRateY);
  PermuYFitted <- zmatrix%*%S%*%PermuY;
  PermuSStot <- sum((PermuY-mean(PermuY))^2);
  PermuSSres <- sum((PermuY-PermuYFitted)^2);
  PermRsquare <- 1-PermuSSres/PermuSStot;
  if(PermRsquare > RSquare) m <- m+1;
 }
 PValue <- m/1000;
 Result <- list(Beta[[lambdaIndex]],Bestlambda,LSerr[lambdaIndex],Penaltyerr[lambdaIndex],PValue,zmatrix,RSquare,covres);
 return(Result)
}