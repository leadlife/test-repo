TangentCheck2 <- function(Xmatrix,tempx,Beta)
{
 if(abs(tempx)<=1)
 { 
  yhat <- Xmatr%*%Beta;
 }
 else
 {
  Funcslope <- function(x){cbind(hermite(x,0),hermite(x,1),hermite(x,2),hermite(x,3),hermite(x,4),hermite(x,1),hermite(x,2),hermite(x,3),hermite(x,4))%*%Beta};
  slope <- grad(Funcslope,sign(tempx));
  tanpointx <- sign(tempx);
  tanpointy <- Funcslope(sign(tempx));
  yhat <- tanpointy+slope*(tempx-tanpointx);
 }
return(yhat);
}