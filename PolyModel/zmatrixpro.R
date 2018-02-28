zmatrixpro <- function(Ndata)
{
 Nsample <- length(Ndata);
 zmatrix <- matrix(0,Nsample,5);
 for(i in 1:Nsample)
 {
  temp <- Ndata[i];
  if(abs(temp)<=1.5)
  {
  zmatrix[i,] <- c(hermite(temp,0),hermite(temp,1),hermite(temp,2),hermite(temp,3),hermite(temp,4));
  }
  else
  {
   zmatrix[i,] <- hermitede(temp);
  }
 }
 return(zmatrix);
}