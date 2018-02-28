DataPro <- function(Origndata,Index,NStart,NRateY)
{
 if(FIndex$Rate[Index] == "Y") 
 {
  Prodata <- log(rev(data.frame(Origndata)[[1]])[c(NStart:(NStart+NRateY-1))]/100+1);
 }
 if(FIndex$Rate[Index] == "B3M") 
 {
  tempdata <- rev(data.frame(Origndata)[[1]])[c(NStart:(NStart+NRateY))]/100;
  Prodata <- log(1+diff(tempdata)*0.25+30/365*tempdata[-c(1)]);
  Prodata <- diff(tempdata)*0.25+30/365*tempdata[-c(1)];
 }
 if(FIndex$Rate[Index] == "B10Y") 
 {
  tempdata <- rev(data.frame(Origndata)[[1]])[c(NStart:(NStart+NRateY))]/100;
  ##Prodata <- log(1+diff(tempdata)*8+30/365*tempdata[-c(1)]);
  Prodata <- diff(tempdata)*8+30/365*tempdata[-c(1)];
 }                                                                                                                                   
 if(FIndex$Rate[Index] == "C") 
 {
  tempdata <- rev(data.frame(Origndata)[[1]])[c(NStart:(NStart+NRateY))]/10000;
  ##Prodata <- log(1+diff(tempdata)*4.5+30/365*tempdata[-c(1)]);
  Prodata <- diff(tempdata)*4.5+30/365*tempdata[-c(1)];
 } 
 if(FIndex$Rate[Index] == "N")
 {
  Prodata <- -diff(log(rev(data.frame(Origndata)[[1]])[c(NStart:(NStart+NRateY))]));
 }
return(Prodata)
}