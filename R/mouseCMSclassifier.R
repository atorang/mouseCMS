#' mouseCMSclassifier
#'
#' \code{mouseCMSclassifier} is a CMS (Consensus Modelcular Subtypes in colorectal cancer)
#'  classifier for mouse organoids.
#'
#' @param data A numeric matrix or data.frame with expression levels. Rows are genes and columns are samples.
#' Ensembl IDs are used.
#' @param perform.log2 If \code{perform.log2 = TRUE}, the data will be log2-transformed.
#' It is a logical determining if data needs log2-transformation or if the data is already log2-transformed.
#' Default is \code{perform.log2 = FALSE}.
#'
#' @return An object of class "\code{data.frame}" containing the classification results.
#' Columns CMS1-4 show the chance of being each subtype, based on an elastic-net model.
#' nearestCMS indicates subtype with highest chance, regardless of the cutoff.
#' predictedCMS, shows subtype of samples that confidently predicted (cuttoff = 0.5).
#' @export
#'
#' @examples
#' #classify a predefined set
#' data("exprs.test")
#' mouseCMSclassifier(data=exprs.test, perform.log2=FALSE)
#'
#'
mouseCMSclassifier <- function(data,perform.log2=FALSE){
  requireNamespace("glmnet", quietly = TRUE)
  data=data.frame(data)
  if(sum(is.na(data))>0){stop("There is missing value (NA) in data")}
  for (i in 1:ncol(data)) {data[,i]=as.numeric(data[,i])}
  if(perform.log2){ data=log2(data+1) }else if(max(data)>30){
    data=log2(data+1)
    warning("Log2-transformation applied on data because values > 30 is observed")
  }
  data.genes=intersect(rownames(exprs.model),rownames(data))
  if(length(data.genes)<0.8*nrow(exprs.model)){
    stop("There is no enough overlapping genes,
         or rownames of data are not Ensembl IDs (like ENSMUSG00000064842)")
  }
  if(length(data.genes)<nrow(exprs.model)){
    warning(paste0("There are some missing genes,
            the classification is less reliable with missing genes:","\n",
                   setdiff(rownames(exprs.model),data.genes)))
  }
    model.test=model
  data=data[rownames(exprs.model),]
  rownames(data)=rownames(exprs.model)
  data[is.na(data)]=0
  for (i in 1:ncol(data)) {
    comb=preprocessCore::normalize.quantiles(as.matrix(data.frame(exprs.model,data[,i])))
    data[,i]=comb[,ncol(comb)]
  }
  prediction=stats::predict(model.test, newx=t(data), interval ="prediction",type="response")
  cls=stats::predict(model.test, newx=t(data), interval ="prediction",type="class")
  Max=apply(prediction, 1, max)
  Predicted=cls
  Predicted[Max<0.5]<-NA
  results=data.frame(prediction,NearestCMS=cls,CMS=Predicted)
  colnames(results)=c("CMS2","CMS3","CMS4","nearestCMS","predictedCMS")
  return(results)
}

