#' Predict proportions using a trained CelMod object
#'
#' This function builds the paired bulk-single cell proportion model from input data, and returns the model object for prediction on new bulk data.
#'
#' @param bdat Matrix of bulk data, with features (genes, proteins, metabolites) as rows and samples as columns.
#' @param classprops Matrix of cell type proportions, with each cell type as a row, and each sample as column. The sample ordering must match those in the bdat bulk matrix.
#' @param numgenevec A vector with range of the number of best-fit features to optimize the model over. Usually 10-20 features is sufficient to model ~10 different cell types
#` @param crossval_times Number of cross-validations to run. For example, crossval_times=5 equates to five-fold cross-validation
#` @param method_type One of "pearson", "spearman", "quantile", or "elasticnet". This parameter determines which fitting model/error function to use. "pearson" uses the Pearson's R to determine the optimal number of features to use, "spearman" uses the Spearman's R, and "quantile" uses a quantile value to identify the max error. "elasticnet" uses an Elastic Net approach for modeling.
#` @param quantileval A parameter to set the quantile value of the max error, if using the "quantile" method
#` @param alphavals A vector of values for the alpha parameter for Elastic Net, over which the model will optimize using cross-validation
#' @return A list containing the following arrays: "cv_preds" with the cross-validation predictions for each value of the numgenevec or alphavals parameter, "cv_errs" with the cross-validation errors for each value of the numgenevec or alphavals parameter, "numgenevec" with the range of best-fit features tested, "cv_bestgenes" with the top features selected for each value of the numgenevec or alphavals parameter, "model" for the trained model, and "modelgenerank" for the final ranking of features used in the model trained on the full data
#' @export
predict_estimates=function(model_list,bdat,typval="pearson") {
  traininds=1:ncol(bdat)
  numgenevec=model_list[["numgenevec"]]
  outprop=list()
  outprop[["numgenevec"]]=numgenevec
  outprop[["proportions"]]=list()
  if (typval=="elasticnet") {
    bdat2=t(bdat)
    for (alphaval in paste0("a_",numgenevec)) {
      y=predict(model_list$model[[alphaval]],bdat2,type="response")
      y2=array(0,dim=dim(y)[c(2,1,3)])
      for (tt in 1:dim(y)[3]) {
        y2[,,tt]=t(y[,,tt])
      }
      outprop[["proportions"]][[alphaval]]=y2
    }
  } else {
    outpred=list()
    for (typ in 1:nrow(classprops)) {
      predval=sweep(sweep(bdat[,traininds,drop=F],1,model_list$model[[typ]][1,],"-"),1,model_list$model[[typ]][2,],"/")
      predval[predval<0]=0
      predval[predval>1]=1
      outpred[[typ]]=predval
    }
    for (keepgen in 1:length(numgenevec)) {
      predmat=matrix(0,nrow=nrow(classprops),ncol=length(traininds))
      for (typ in 1:nrow(classprops)) {
        predmat[typ,]=apply(outpred[[typ]][model_list$modelgenerank[1:numgenevec[keepgen],typ],],2,function(x){return(mean(x,na.rm=T))})
      }
      rownames(predmat)=rownames(model_list$cv_preds[[1]])
      outprop[["proportions"]][[keepgen]]=sweep(predmat,2,colSums(predmat),"/")
    }
  }
  return(outprop)
}
