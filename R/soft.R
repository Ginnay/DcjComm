
soft = function( x, A ){
  if (sum( abs(cbind(A)) )==0){
    y <- x
  }else{
    y<- (abs(x) - A)
    y[y < 0] <- 0
    y <- sign(x)*y
  }
  return (y)
   
}


