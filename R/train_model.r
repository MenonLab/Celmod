#' Train a Celmod model
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
#' @export#' @export

train_model=function(bdat,classprops,numgenevec=3:10,crossval_times=5,seedval=1,method_type="pearson",quantileval=0.9,alphavals=seq(from=0,to=1,length.out=11)) {
  ###select training set for CV###
  keepcells=rep(1:crossval_times,length.out=ncol(bdat))
  set.seed(seedval)
  keepcells=sample(keepcells)
  cv_preds=list()
  if (method_type=="elasticnet") {
    bdat=t(bdat)
    allcombs=expand.grid(alphavals,1:100)
    allcombs=paste0(allcombs[,1],"_",allcombs[,2])
    for (ii in allcombs) {
      cv_preds[[ii]]=matrix(NA,ncol=nrow(classprops),nrow=ncol(classprops))
    }
    for (alphaval in alphavals) {
      print(paste0("running cross-validation for alpha ",alphaval))
      for (cvalind in 1:crossval_times) {
        traininds=which(keepcells!=cvalind)
        testinds=which(keepcells==cvalind)
        y=glmnet::glmnet(bdat[traininds,],y=t(classprops[,traininds]),family="multinomial",alpha=alphaval)
        out1=predict(y,bdat[testinds,,drop=F],type="response")
        for (jj in 1:dim(out1)[3]) {
          cv_preds[[paste0(alphaval,"_",jj)]][testinds,]=(out1[,,jj])      
        } 
      }
    }
    ###calculate correlations for each parameter set###
    cv_errs=matrix(NA,nrow=nrow(classprops),ncol=length(cv_preds))
    rownames(cv_errs)=rownames(classprops)
    colnames(cv_errs)=names(cv_preds)
    for (ii in 1:ncol(cv_errs)) {
      cv_errs[,ii]=diag(cor(cv_preds[[ii]],t(classprops)))
    }
    ddd=apply(cv_errs,2,function(x){return(min(x,na.rm=T))})
    numgene=which.max(apply(cv_errs,2,function(x){return(min(x,na.rm=T))})) ###max of minimum value
    
    allmods=list()
    traininds=1:nrow(bdat)
    for (alphaval in alphavals) {
        y=glmnet(bdat[traininds,],y=t(classprops[,traininds]),family="multinomial",alpha=alphaval)
        allmods[[paste0("a_",alphaval)]]=y
    }
    outlist=list()
    outlist[["cv_preds"]]=cv_preds
    outlist[["cv_errs"]]=cv_errs
    outlist[["numgenevec"]]=alphavals
    outlist[["cv_bestgenes"]]=numgene
    outlist[["model"]]=allmods
    outlist[["modelgenerank"]]=NULL 
  } else {
    for (ii in 1:length(numgenevec)) {
      cv_preds[[ii]]=matrix(NA,nrow=nrow(classprops),ncol=ncol(bdat))
    }
    for (cvalind in 1:crossval_times) {
      print(paste0("running cross-validation ",cvalind))
      traininds=which(keepcells!=cvalind)
      testinds=which(keepcells==cvalind)
      allmods=list()
      trainpred=list()
      testpred=list()
      for (typ in 1:nrow(classprops)) {
        y=apply(bdat[,traininds],1,function(x){outval=summary(lm(x~classprops[typ,traininds]))$coefficients;return(c(outval[1,1],outval[2,1]))})
        allmods[[typ]]=y
        predval=sweep(sweep(bdat[,traininds,drop=F],1,allmods[[typ]][1,],"-"),1,allmods[[typ]][2,],"/")
        predval[predval<0]=0
        predval[predval>1]=1
        trainpred[[typ]]=predval
        predval=sweep(sweep(bdat[,testinds,drop=F],1,allmods[[typ]][1,],"-"),1,allmods[[typ]][2,],"/")
        predval[predval<0]=0
        predval[predval>1]=1
        testpred[[typ]]=predval
      }
      ###order genes by fit for each type
      traingene_rank=c()
      for (typ in 1:nrow(classprops)) {
        if (method_type=="pearson") {
          genevals=cor(t(trainpred[[typ]]),(classprops[typ,traininds]))
          traingene_rank=cbind(traingene_rank,order(-abs(genevals)))
        }
        if (method_type=="spearman") {
          genevals=cor(t(trainpred[[typ]]),(classprops[typ,traininds]),method="spearman")
          traingene_rank=cbind(traingene_rank,order(-abs(genevals)))
        }
        if (method_type=="quantile") {
          genevals=abs(sweep(trainpred[[typ]],2,classprops[typ,traininds],"-"))
          genevals=1-apply(genevals,1,function(x){quantile(x,quantileval,na.rm = T)})
          traingene_rank=cbind(traingene_rank,order(-abs(genevals)))
        } 
      }
      ####apply ranked gene values to test list###
      for (keepgen in 1:length(numgenevec)) {
        predmat=matrix(0,nrow=nrow(classprops),ncol=length(testinds))
        for (typ in 1:nrow(classprops)) {
          predmat[typ,]=apply(testpred[[typ]][traingene_rank[1:numgenevec[keepgen],typ],,drop=F],2,function(x){return(mean(x,na.rm=T))})
        }
        predmat=sweep(predmat,2,colSums(predmat),"/")
        
        cv_preds[[keepgen]][,testinds]=predmat
        rownames(cv_preds[[keepgen]])=rownames(classprops)
      }
    }
    ###err values###
    cv_errs=matrix(NA,nrow=nrow(classprops),ncol=length(numgenevec))
    for (keepgen in 1:length(numgenevec)) {
      if (method_type=="pearson") {
        cor1=cor(t(classprops),t(cv_preds[[keepgen]]),method="pearson")
        cv_errs[,keepgen]=diag(cor1)
      }
      if (method_type=="spearman") {
        cor1=cor(t(classprops),t(cv_preds[[keepgen]]),method="spearman")
        cv_errs[,keepgen]=diag(cor1)
      }
      if (method_type=="quantile") {
        cor1=1-abs(cv_preds[[keepgen]]-classprops)
        cv_errs[,keepgen]=apply(cor1,1,function(x){quantile(x,quantileval,na.rm=T)})
      }
    }
    numgene=which.max(apply(cv_errs,2,function(x){return(min(x,na.rm=T))}))
    
    ###fullmodel###
    trainpred=list()
    allmods=list()
    traininds=1:ncol(bdat)
    for (typ in 1:nrow(classprops)) {
      y=apply(bdat[,traininds],1,function(x){outval=summary(lm(x~classprops[typ,traininds]))$coefficients;return(c(outval[1,1],outval[2,1]))})
      allmods[[typ]]=y
      predval=sweep(sweep(bdat[,traininds,drop=F],1,allmods[[typ]][1,],"-"),1,allmods[[typ]][2,],"/")
      predval[predval<0]=0
      predval[predval>1]=1
      trainpred[[typ]]=predval
    }
    ###order genes by fit for each type
    traingene_rank=c()
    for (typ in 1:nrow(classprops)) {
      if (method_type=="pearson") {
        genevals=cor(t(trainpred[[typ]]),classprops[typ,traininds])
        traingene_rank=cbind(traingene_rank,order(-abs(genevals)))
      }
      if (method_type=="spearman") {
        genevals=cor(t(trainpred[[typ]]),classprops[typ,traininds],method="spearman")
        traingene_rank=cbind(traingene_rank,order(-abs(genevals)))
      }
      if (method_type=="quantile") {
        genevals=abs(sweep(trainpred[[typ]],2,classprops[typ,traininds],"-"))
        genevals=1-apply(genevals,1,function(x){quantile(x,0.9,na.rm=T)})
        traingene_rank=cbind(traingene_rank,order(-abs(genevals)))
      } 
    }
    outlist=list()
    outlist[["cv_preds"]]=cv_preds
    outlist[["cv_errs"]]=cv_errs
    outlist[["numgenevec"]]=numgenevec
    outlist[["cv_bestgenes"]]=numgene
    outlist[["model"]]=allmods
    outlist[["modelgenerank"]]=traingene_rank
  }
  return(outlist)
}