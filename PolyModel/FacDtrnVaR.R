FacDtrnVaR <- function(ndata,Facmin,Index,NRateY)
{
##input
 ##ndata:  price(order: newest to oldest)
 ##Facmin:  threshold(Tail beinning point)
##output
 ##nVaR: VaR(@ 99%(pareto Tail Distribution),84%,50%,16%,1%(Normal Distribution));
Return <- DataPro(ndata,Index,NRateY);
average <- mean(Return);
stdv <- sd(Return);
NReturn <- -(Return-mean(Return))/sd(Return);
TailSample <- subset(NReturn,NReturn>Facmin);
LenSample <- length(NReturn);
LenTail <- length(TailSample);
TailProb <- 1-(((LenSample-LenTail+1):LenSample)-0.5)/LenSample;
OrderTail <- sort(TailSample);
OrderTailProb <- cbind(OrderTail,TailProb);
OrderLogTailProb <- log(OrderTailProb);
OrderLogTailProb <- data.frame(OrderLogTailProb);
Tailfit <- lm(TailProb~OrderTail,OrderLogTailProb);
##x11();
##plot(OrderLogTailProb$OrderTail,OrderLogTailProb$TailProb);
##lines(OrderLogTailProb$OrderTail,Tailfit$fitted,col=4,lwd=2);

nConstant <- exp(summary(Tailfit)$coefficients[1]); 
nPowerLaw <- -summary(Tailfit)$coefficients[2];
##nConstant 
##nPowerLaw 

xx <- seq(Facmin,max(TailSample),0.2);
##x11();
##plot(xx,nConstant*((xx)^(-nPowerLaw)));
##plot(OrderTailProb);
##lines(xx,nConstant*((xx)^(-nPowerLaw)),col="red",lwd=2);
VaR99 <- -(nConstant/(1-0.99))^(1/nPowerLaw);

VaR84 <- average+stdv*qnorm(1-0.84,0,1);
VaR50 <- average+stdv*qnorm(1-0.50,0,1);
VaR16 <- average+stdv*qnorm(1-0.16,0,1);
VaR1 <- average+stdv*qnorm(1-0.1,0,1);
VaR <- rbind(VaR99,VaR84,VaR50,VaR16,VaR1);
return(VaR)
}