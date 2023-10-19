#' Normalize the input data
#' @param params  The input params.
#' @export
#'
Normalize <- function (params){
U <- params$U
V <- params$V
NormV<-params$NormV
Norm<-params$Norm
nSmp <- dim(V)[1]
mFea <- dim(U)[1]

if (Norm == 2){
  if (NormV){
    norms <- sqrt(apply(V*V,2,sum))
    norms <- as.matrix(norms)
    V <- V/(matrix(t(norms),nrow=nSmp,ncol=length(norms),byrow=T))
    U <- U*(matrix(t(norms),nrow=mFea,ncol=length(norms),byrow=T))
  } else{
    norms <- sqrt(apply(U*U,2,sum))
    norms <- as.matrix(norms)
    U <- U/(matrix(t(norms),nrow=mFea,ncol=length(norms),byrow=T))
    V <- V*(matrix(t(norms),nrow=nSmp,ncol=length(norms),byrow=T))
  }
}else
  if (NormV){
    norms <- apply(abs(V),2,sum)
    norms <- as.matrix(norms)
    V <- V/(matrix(t(norms),nrow=nSmp,ncol=length(norms),byrow=T))
    U <- U*(matrix(t(norms),nrow=mFea,ncol=length(norms),byrow=T))
  } else{
    norms <- apply(abs(U),2,sum)
    norms <- as.matrix(norms)
    U <- U/(matrix(t(norms),nrow=mFea,ncol=length(norms),byrow=T))
    V <- V*(matrix(t(norms),nrow=nSmp,ncol=length(norms),byrow=T))
  }

 Normalize_Results <- list(U=U, V=V)
 #Normalize_Results <- list(NormalizeUV_Results = NormalizeUV_Results)

 return(Normalize_Results)
 #return(V)
}
