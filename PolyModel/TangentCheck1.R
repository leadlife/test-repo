TangentCheck1 <- function(tempx,Beta)
{
  Funcslope <- function(x){cbind(hermite(x,0),hermite(x,1),hermite(x,2),hermite(x,3),hermite(x,4))%*%Beta};
  slope <- grad(Funcslope,tempx);
  tanpointx <- tempx;
  tanpointy <- Funcslope(tempx);
  yhat <- tanpointy+slope*(tempx-tanpointx);
  intercept <- tanpointy-slope*tanpointx;
  result <- list(yhat,slope,intercept);
  return(result);
}