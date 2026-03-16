library(lme4)
library(randomForest)
library(e1071)
library(glmertree)
library(dplyr)
library(tibble)
library(tidyr)
library(stringr)
library(tictoc)
library(writexl)

tree=function(x4,x5,x6){
  case_when(
    x4<=1 ~ case_when(
      abs(x6)>=3 ~ case_when(
        x5<=5 ~ 20,
        x5>5  ~ -2 
      ),
      abs(x6)<3 ~ 2
    ),
    x4>1 ~ case_when(
      x5 >= -1 ~ 0,
      x5 < -1  ~ case_when(
        x6<=0 ~ 1,
        x6>0  ~ -1
      )
    )
  )
}

fixed=function(x1,x2,x3,x4,x5,x6,x7,alpha,beta){
  alpha*(x1*2-3*x2-x2*x3^2)+beta*tree(x4,x5,x6)
}

sim_gmexgb=function(N,n_train,n_test,alpha,beta,gama2,delta2){
  n=n_train+n_test
  group=as.factor(rep(1:N,each=n))
  train_test=rep(rep(c("train","test"),c(n_train,n_test)),times=N)
  
  x1=runif(N*n,-1,1)
  x2=runif(N*n,-1,1)
  x3=rweibull(N*n,3)
  x4=runif(N*n,-3,3)
  x5=runif(N*n,-6,6)
  x6=runif(N*n,-5,5)
  x7=runif(N*n,-4,4)
  
  fix=fixed(x1,x2,x3,x4,x5,x6,x7,alpha,beta)
  
  b0=rep(rnorm(N,0,sd=sqrt(gama2)),each=n)
  b1=rep(rnorm(N,0,sd=sqrt(delta2)),each=n)
  rand=b0+b1*x1
  
  eta=fix+rand
  mu=plogis(eta)
  y=as.factor(rbinom(N*n,size=1,prob=mu))
  
  out=data.frame(
    group,train_test,x1,x2,x3,x4,x5,x6,x7,y,b0,b1,fix,rand,eta,mu
  )
  out
}

sim_pref=function(N,n_train,n_test,alpha,beta,gama2,delta2){
  pmad=NA
  pmcr=NA
  
  s1=sim_gmexgb(N,n_train,n_test,alpha,beta,gama2,delta2)
  
  trainset=s1 |> 
    filter(train_test=="train")
  testset=s1 |> 
    filter(train_test=="test")
  
  #######################glmm xgb
  if(delta2 == 0){
    fit_gmexgb=try(
      gmexgb(y=trainset$y, cov=trainset,
             group = trainset$group, 
             xnam=paste0("x",1:7),
             family='binomial'),
      silent = T
    )
  } else {
    fit_gmexgb=try(
      gmexgb(y=trainset$y, cov=trainset,
             group = trainset$group, 
             xnam=paste0("x",1:7),
             znam="x1",
             family='binomial'),
      silent = T
    )
  }
  
  # fit_gmexgb
  if(! inherits(fit_gmexgb, "try-error")){
    pred_resp=predict.gmexgb(fit_gmexgb,testset,group=testset$group,type='response')
    # table(pred_resp)
    yhat=pred_resp |> as.integer() |> as.factor()
    # table(yhat)
    y=testset$y
    # table(y)
    # table(yhat,y)
    
    pmcr=mean(yhat != y)
    # pmcr
    
    ############# roc
    pred_mu=predict.gmexgb(fit_gmexgb,testset,group=testset$group,type='mu')
    mu=testset$mu
    # plot(mu,pred_mu)
    pmad=mean(abs(mu-pred_mu))
    # pmad
  }
  out=c(pmad=pmad,pmcr=pmcr)
}

sim_pref_total=function(N,n_train,n_test,alpha,beta,gama2,delta2,
                        eta=0.3,nrounds=100,max_depth=6,subsample=1,lambda=1,alpha_xgb=0){
  
  s1=sim_gmexgb(N,n_train,n_test,alpha,beta,gama2,delta2)
  
  trainset=s1 |> 
    filter(train_test=="train")
  testset=s1 |> 
    filter(train_test=="test")
  
  y=testset$y
  # table(y)
  mu=testset$mu
  
  ########################## xgb
  predictors = trainset |> 
    select(group,x1:x7) |> 
    data.matrix()
  label = as.factor(as.numeric(trainset$y)-1)
  fit_xgb = xgboost(
    data=predictors, label=label,
    params=list(
      subsample=subsample, eta=eta, max_depth=max_depth,
      lambda=lambda, alpha=alpha_xgb
      # eta=0.3,max_depth=6,subsample=1,lambda=1,alpha=0
    ),
    nrounds=nrounds,
    verbose = 0, 
    objective="binary:logistic"
  )
  # fit_xgb
  newdata=testset |> 
    select(group,x1:x7) |> 
    data.matrix()
  pred_resp=predict(fit_xgb,newdata)
  # table(pred_resp)
  yhat=(pred_resp >= 0.5) |> as.integer() |> as.factor()
  # table(yhat)
  # table(yhat,y)
  
  pmcr_xgb=mean(yhat != y)
  # pmcr_xgb
  
  ############# roc
  pred_mu=pred_resp
  # plot(mu,pred_mu)
  pmad_xgb=mean(abs(mu-pred_mu))
  # pmad_xgb
  
  #######################gmexgb
  # fit_gmexgb
  if(delta2 == 0){
    fit_gmexgb=try(
      gmexgb(
        y=trainset$y, cov=trainset,
        group = trainset$group, 
        xnam=paste0("x",1:7),
        family='binomial',
        eta=eta,nrounds=nrounds,
        max_depth=max_depth,subsample=subsample,
        lambda=lambda,alpha=alpha_xgb
      ),
      silent = T
    )
  } else {
    fit_gmexgb=try(
      gmexgb(
        y=trainset$y, cov=trainset,
        group = trainset$group, 
        xnam=paste0("x",1:7),
        znam="x1",
        family='binomial',
        eta=eta,nrounds=nrounds,
        max_depth=max_depth,subsample=subsample,
        lambda=lambda,alpha=alpha_xgb
      ),
      silent = T
    )
  }
  
  if(! inherits(fit_gmexgb, "try-error")){
    pred_resp=predict.gmexgb(fit_gmexgb,testset,group=testset$group,type='response')
    # table(pred_resp)
    yhat=pred_resp |> as.integer() |> as.factor()
    # table(yhat)
    # table(yhat,y)
    
    pmcr_gmexgb=mean(yhat != y)
    # pmcr_gmexgb
    
    ############# roc
    pred_mu=predict.gmexgb(fit_gmexgb,testset,group=testset$group,type='mu')
    # plot(mu,pred_mu)
    pmad_gmexgb=mean(abs(mu-pred_mu))
    # pmad_gmexgb
  } else {
    pmcr_gmexgb=pmad_gmexgb=NA
  }
  
  ########################## glm
  fit_glm=glm(
    y~group+x1+x2+x3+x4+x5+x6+x7,
    family = "binomial",
    data = trainset
  )
  # summary(fit_glm)
  pred_resp=predict(fit_glm,testset,type='response')
  # table(pred_resp)
  yhat=(pred_resp >= 0.5) |> as.integer() |> as.factor()
  # table(yhat)
  # table(yhat,y)
  
  pmcr_glm=mean(yhat != y)
  # pmcr_glm
  
  ############# roc
  pred_mu=pred_resp
  # plot(mu,pred_mu)
  pmad_glm=mean(abs(mu-pred_mu))
  # pmad_glm
  
  ########################## glmer
  if(delta2 == 0){
    fit_glmer=glmer(
      y~(1|group)+x1+x2+x3+x4+x5+x6+x7,
      family = "binomial",
      data = trainset
    )
  } else {
    fit_glmer=glmer(
      y~(1|group)+(x1|group)+x2+x3+x4+x5+x6+x7,
      family = "binomial",
      data = trainset
    )
  }
  # summary(fit_glmer)
  pred_resp=predict(fit_glmer,testset,type='response')
  # table(pred_resp)
  yhat=(pred_resp >= 0.5) |> as.integer() |> as.factor()
  # table(yhat)
  # table(yhat,y)
  
  pmcr_glmer=mean(yhat != y)
  # pmcr_glmer
  
  ############# roc
  pred_mu=pred_resp
  # plot(mu,pred_mu)
  pmad_glmer=mean(abs(mu-pred_mu))
  # pmad_glmer
  
  ########################## rf
  fit_rf=randomForest(
    y~group+x1+x2+x3+x4+x5+x6+x7,
    data = trainset
  )
  # fit_rf
  pred_resp=predict(fit_rf,testset,type='response')
  # table(pred_resp)
  yhat=pred_resp
  # table(yhat)
  # table(yhat,y)
  
  pmcr_rf=mean(yhat != y)
  # pmcr_rf
  
  ############# roc
  pred_mu=predict(fit_rf,testset,type='prob')[,"1"]
  # plot(mu,pred_mu)
  pmad_rf=mean(abs(mu-pred_mu))
  # pmad_rf
  
  ########################## svm
  fit_svm=svm(
    y~group+x1+x2+x3+x4+x5+x6+x7,
    data = trainset,
    probability=T
  )
  # fit_svm
  pred_resp=predict(fit_svm,testset)
  # table(pred_resp)
  yhat=pred_resp
  # table(yhat)
  # table(yhat,y)
  
  pmcr_svm=mean(yhat != y)
  # pmcr_svm
  
  ############# roc
  pred=predict(fit_svm,testset,probability=T)
  pred_mu=attr(pred,"probabilities")[,"1"]
  # plot(mu,pred_mu)
  pmad_svm=mean(abs(mu-pred_mu))
  # pmad_svm
  
  ########################## glmertree
  if(delta2 == 0){
    fit_glmertree=glmertree(
      y~1 | group | x1+x2+x3+x4+x5+x6+x7,
      family = "binomial",
      data = trainset
    )
  } else {
    fit_glmertree=glmertree(
      y~x1 | group | x2+x3+x4+x5+x6+x7,
      family = "binomial",
      data = trainset
    )
  }
  # fit_glmertree
  # plot(fit_glmertree, which = "all") # default behavior, may also be "tree" or "ranef" 
  # coef(fit_glmertree)
  # ranef(fit_glmertree)
  pred_resp=predict(fit_glmertree,testset, type = "response")
  # table(pred_resp)
  yhat=(pred_resp >= 0.5) |> as.integer() |> as.factor()
  # table(yhat)
  # table(yhat,y)
  
  pmcr_glmertree=mean(yhat != y)
  # pmcr_glmertree
  
  ############# roc
  pred_mu=pred_resp
  # plot(mu,pred_mu)
  pmad_glmertree=mean(abs(mu-pred_mu))
  # pmad_glmertree
  
  ####################### gmerf
  if(delta2 == 0){
    fit_gmerf=try(
      gmerf(y=trainset$y, cov=trainset,
            group = trainset$group, 
            xnam=paste0("x",1:7),
            family='binomial'),
      silent = T
    )
  } else {
    fit_gmerf=try(
      gmerf(y=trainset$y, cov=trainset,
            group = trainset$group, 
            xnam=paste0("x",1:7),
            znam="x1",
            family='binomial'),
      silent = T
    )
  }
  
  if(! inherits(fit_gmerf, "try-error")){
    pred_resp=predict.gmerf(fit_gmerf,testset,group=testset$group,type='response')
    # table(pred_resp)
    yhat=pred_resp |> as.integer() |> as.factor()
    # table(yhat)
    # table(yhat,y)
    
    pmcr_gmerf=mean(yhat != y)
    # pmcr_gmerf
    
    ############# roc
    pred_mu=predict.gmerf(fit_gmerf,testset,group=testset$group,type='mu')
    # plot(mu,pred_mu)
    pmad_gmerf=mean(abs(mu-pred_mu))
    # pmad_gmerf
  } else {
    pmcr_gmerf=pmad_gmerf=NA
  }
  
  out=c(
    pmad_glm=pmad_glm,pmcr_glm=pmcr_glm,
    pmad_glmer=pmad_glmer,pmcr_glmer=pmcr_glmer,
    pmad_rf=pmad_rf,pmcr_rf=pmcr_rf,
    pmad_svm=pmad_svm,pmcr_svm=pmcr_svm,
    pmad_glmertree=pmad_glmertree,pmcr_glmertree=pmcr_glmertree,
    pmad_gmerf=pmad_gmerf,pmcr_gmerf=pmcr_gmerf,
    pmad_xgb=pmad_xgb,pmcr_xgb=pmcr_xgb,
    pmad_gmexgb=pmad_gmexgb,pmcr_gmexgb=pmcr_gmexgb
  )
}

table_summary=function(tb){
  m1=colMeans(tb,na.rm = T)
  v1=matrixStats::colVars(tb,na.rm = T)
  
  d1=tibble(stat="Mean",id=names(m1),value=m1)
  d1
  d2=tibble(stat="Var",id=names(m1),value=v1)
  d2
  
  d3=bind_rows(d1,d2)
  d3
  
  d4=d3 |> 
    separate_wider_delim(id,delim="_",names=c("measure","Algorithm")) |> 
    mutate(measure=str_to_upper(measure))
  d4
  
  d5=d4 |> 
    pivot_wider(
      names_from = c(stat,measure),
      names_sep = "_",
      values_from = value
    ) |> 
    relocate(Algorithm, Mean_PMAD, Var_PMAD, Mean_PMCR, Var_PMCR)
  d5
}


