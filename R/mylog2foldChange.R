#' calculate foldchange for each celltype
#' @param inData a dataframe of gene expression
#' @param cell.type the cell type which you want to calculate the foldchange
#' @param method default is "median",eg "median" or "mean"
#' @param probs the quantile for median calculate
#' @return the foldchange value of each celltype in \code{cellwave object}
#' @importFrom stats median

mylog2foldChange.diy<-function(inData, cell.type, method="median", probs = 0.75) # median or mean
{
  re.list <- list()
  cell.type <- cell.type
  print(cell.type)
  for (c in cell.type) {
    print(c)
    cell_fc_labels<-colnames(inData)
    cell_fc_labels[cell_fc_labels==c]<-0
    cell_fc_labels[cell_fc_labels!='0']<-1

    classLabel <- cell_fc_labels
    sampleIdsCase<-which(classLabel==0);#0 tumer
    sampleIdsControl<-which(classLabel==1);#1 normal
    probeFC <- as.numeric(apply(inData, 1, function(x){
      # print(x)
      if(method=="mean"){
        (mean(as.numeric(x[sampleIdsCase]))+1)/(mean(as.numeric(x[sampleIdsControl]))+1);
      }else if(method=="median"){
        quantile(as.numeric(x[sampleIdsCase]), probs = probs, names=FALSE)

        a.tmp <- quantile(as.numeric(x[sampleIdsCase]), probs = probs, names=FALSE)+1
        b.tmp <- quantile(as.numeric(x[sampleIdsControl]), probs = probs, names=FALSE)+1
        a.tmp/b.tmp
      }
    }))

    probeFC<-log(probeFC,base=2);
    fc_res<-probeFC;

    fc_res[is.infinite(fc_res)]<-0
    fc_res[is.na(fc_res)]<-0
    res = data.frame(gene_id = as.vector(rownames(inData)), log2fc = fc_res)
    # print(dim(res))
    re.list <- c(re.list, list(res))
  }
  names(re.list) <- cell.type
  return(re.list)
}






