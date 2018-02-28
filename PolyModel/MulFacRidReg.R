MulFacRidReg <- function(lambda,Alldata,IndexDep,NFactBegin,NFactEnd,NStart,NRateY)
{
 result <- list();
 LSErr <- rep(0,length(N));
 PenErr <- rep(0,length(N));
 ErrRatio <- rep(0,length(N));
 ErrDiff <- rep(0,length(N));
 Pvalue <- rep(0,length(N));
 RSquared <- rep(0,length(N));
 Covres <- rep(0,length(N));
 FactdataPro <- list();
 DepdataPro <- list();
 for (i in NFactBegin:NFactEnd)
 {
  Depdatatemp <- Alldata[[IndexDep]][index(Alldata[[i]])];
  Depdata <- DataPro(Depdatatemp,IndexDep,NStart,NRateY);
  Facdata <- DataPro(Alldata[[i]],i,NStart,NRateY);
  NDepLogReturn <- Depdata-mean(Depdata);
  NFacLogReturn <- (Facdata-mean(Facdata))/sd(Facdata);
  DepdataPro[[i]] <- NDepLogReturn;
  FactdataPro[[i]] <- NFacLogReturn;
  result[[i]] <- OneFacRidgReg(NDepLogReturn,NFacLogReturn,lambda);
  LSErr[i] <- result[[i]][[3]];
  PenErr[i] <- result[[i]][[4]]; 
  ErrRatio[i] <- abs(result[[i]][[4]]/result[[i]][[3]]-1);
  ErrDiff[i] <- abs(LSErr[i]-PenErr[i]);
  Pvalue[i] <- result[[i]][[5]];
  RSquared[i] <- result[[i]][[7]];
  Covres[i] <- result[[i]][[8]];
  ##print(i);
 }
 VOutput <- list(LSErr,ErrDiff,ErrRatio,DepdataPro,FactdataPro,result,Pvalue,RSquared,Covres);
 return(VOutput);
}