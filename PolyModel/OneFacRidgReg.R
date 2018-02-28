OneFacRidgReg <- function(NDepdata, NFacdata, lambda)
{
 y <- NDepdata;
 x <- NFacdata;
 ##x <- NIXICETFLogReturn;
 ##Penamatrix <- diag(1,5);
 Penamatrix <- matrix(c(0,0,0,0,0, 0,0,0,0,0, 0,0,4,0,0, 0,0,0,9,0, 0,0,0,0,16),5,5);
 Nsample <- length(y);
 LSerr <- rep(0,length(lambda));
 Penaltyerr <- rep(0,length(lambda));
 Differr <- rep(0,length(lambda));
 tempx <- x;
 tempy <- y;
 zmatrix <- zmatrixpro(tempx);
 Beta <- list();   
 for (i in 1:length(lambda))
 {
   S <- matrix.inverse(t(zmatrix)%*%zmatrix+lambda[i]*Penamatrix)%*%t(zmatrix);
   Beta[[i]] <- S%*%tempy;
   Fittedvalue <- zmatrix%*%Beta[[i]];
   LSerr[i] <- sum((Fittedvalue - tempy)^2);
   Penaltyerr[i] <- lambda[i]*diag(Penamatrix)^2%*%Beta[[i]]^2; 
   
 }
 ratio <- abs(Penaltyerr/LSerr-1);
 temp <- which(ratio <= max(head(sort(ratio),10)));
 lambdaIndex <- temp[which.min(LSerr[temp])];
 Bestlambda <- lambda[lambdaIndex];
 ##Bestlambda <- 0;
 SSres <- sum((zmatrix%*%Beta[[lambdaIndex]]- tempy)^2);
 covres <- SSres/(length(y)-5);
 SStot <- sum((y-mean(y))^2);
 RSquare <- 1-SSres/SStot;
 m <- 0;
 S <- matrix.inverse(t(zmatrix)%*%zmatrix+Bestlambda*Penamatrix)%*%t(zmatrix);
##Pvalue
 RsquarePerm <- vector(mode = "numeric", length = 1000)
 for(i in 1:1000)
 {
  PermuY <- sample(tempy,NRateY);
  PermuYFitted <- zmatrix%*%S%*%PermuY;
  PermuSStot <- sum((PermuY-mean(PermuY))^2);
  PermuSSres <- sum((PermuY-PermuYFitted)^2);
  RsquarePerm[i] <- 1-PermuSSres/PermuSStot;   
 }
    
 ##Rmean <- mean(RsquarePerm);
 ##Rsd <- sd(RsquarePerm);
 R <- 1/(1-RsquarePerm);
 Rthreshold <- sort(R, decreasing = T)[100];
   
 Rtail <- sort(R, decreasing = T)[5:100];
 
 nRtail <- length(Rtail);
 tailprob <- 1 - ((1000-nRtail+1):1000-0.5)/1000;
 Ordertail <- sort(Rtail)
 Ordertailprob <- cbind(Ordertail, tailprob);
 Orderlogtailprob <- log(Ordertailprob);
 Orderlogtailprob <- as.data.frame(Orderlogtailprob)
    
 tailfit <- lm(tailprob ~ Ordertail, data = Orderlogtailprob);
    
    
 nConstant <- as.numeric(exp(tailfit$coefficients[1]));
 nPowerLaw <- as.numeric(-tailfit$coefficients[2]);
    
    
 if (1/(1-RSquare) > Rthreshold)
 {
  PValue <- nConstant * (1/(1-RSquare))^(-nPowerLaw)
 } 
 else 
 {
  PValue <- 1 - pnorm(1/(1-RSquare), mean = mean(R), sd = sd(R))
 }
 #for (j in 1:1000)
 #{
 #PermuY <- sample(tempy,NRateY);
 #PermuYFitted <- zmatrix%*%S%*%PermuY;
 #PermuSStot <- sum((PermuY-mean(PermuY))^2);
 #PermuSSres <- sum((PermuY-PermuYFitted)^2);
 #PermRsquare <- 1-PermuSSres/PermuSStot;
 #if(PermRsquare > RSquare) m <- m+1;
 #}
 #PValue <- m/1000;
 Result <- list(Beta[[lambdaIndex]],Bestlambda,LSerr[lambdaIndex],Penaltyerr[lambdaIndex],PValue,zmatrix,RSquare,covres);
 return(Result)
}