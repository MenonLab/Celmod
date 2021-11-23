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
