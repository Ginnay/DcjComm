fisher_score <- function(Tup,Ltf,genes)
{
  a <- length(intersect(Tup,Ltf))
  b <- length(Tup)-a
  c <- length(Ltf)-a
  d <- length(genes)-(a+b+c)
  
  mat <- matrix(c(a,c,b,d),2,2)
  score <- 1 - fisher.test(mat,alternative="greater")$p.value
  return(score)
}