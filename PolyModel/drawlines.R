drawlines <- function(datax,datay) 
{
 datayy <- matrix(0,length(datay),1);
 for (i in 1:length(datay))
 {
  datayy[i] <- datay[[i]];
 }
 data <- data.frame(datax,datay);
 j <- order(data$datax);
 Result <- lines(data$datax[j],data$datay[j],col="blue",lwd=2);
 return(Result);
}