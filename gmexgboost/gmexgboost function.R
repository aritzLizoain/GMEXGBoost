gmexgb= function (y, cov, group, xnam=NULL, znam=NULL, family='binomial', bizero=NULL, 
                  itmax=30, toll=0.02,
                  eta=0.3,nrounds=100,max_depth=6,subsample=1,lambda=1,alpha=0) {
  #subjects: 
  #-y=vector of responses
  #-cov=data frame with covariates of each statistical unit
  #-group=A factor vector that, for each statistical unit, indicates the group to which it belongs
  #-xnam=vector with names of covariates to be used in xgb
  #-znam=vector with names of covariates to be used in random effects
  #-bizero=A matrix in which each column contains the coefficients of the random effects
  #	   the first value of each column is the intercept, the others are the covariates of znam
  #I assume that group and bizero are consistent, that is? b[,i] matches levels(group)[i]
  
  ######################################
  ####	STEP 1:Initialization    ######
  ######################################
  
  N <- length(y) #number of observations
  n=length(levels(group)) #number of groups
  q <- length(znam)+1	# number of covariates + random intercept
  Zi=NULL
  z.not.null=!(is.null(znam)) #check if there are covariates included in the random effects
  if (z.not.null) Zi=cov[znam] #random effects covariates
  Zi.int= cbind(rep(1,N),Zi) #random intercept + random effects covariates
  
  #Initialize bi to 0 if it is NULL
  if( is.null(bizero) ){
    bi <- NULL
    for(i in 1:n) bi=cbind(bi,rep(0,q))
  }	
  if( !is.null(bizero) ) bi=bizero
  lev=levels(group) #nomi dei gruppi
  bi=data.frame(bi)
  names(bi)=lev
  all.bi=list()  #the b_i of each iteration
  all.bi[[1]]=bi
  
  #If xnam is NULL, assume that all covariates will be used
  if(is.null(xnam)) xnam=names(cov)
  
  #group must be a factor, otherwise it gives an error
  if(!is.factor(group)) stop('Argument "group" must be a factor')
  
  if(z.not.null){
    glmer.formula=as.formula(paste("y ~ ( 1+", paste(znam, collapse= "+"), " | group )"))
  }
  if(!z.not.null){glmer.formula=as.formula(paste("y ~ ( 1 | group)"))}
  xgb.formula=as.formula(paste("target ~ ", paste(xnam, collapse= "+")))
  
  
  ##################################################
  ####	STEP 2: GLM to initialize mu_ij    #######
  ##################################################
  
  glm.formula=as.formula(paste("y ~ ", paste(xnam, collapse= "+")))
  glm.data= cbind(y,cov[xnam])
  glm.0=glm(glm.formula, data=glm.data, family=family)
  mu.ij.0=glm.0$fitted.values
  eta.est=glm.0$family$linkfun(mu.ij.0) #glm estimates
  linkf=glm.0$family$linkfun
  linkinv=glm.0$family$linkinv
  
  
  ####################################################
  ####	STEP 3-4: Iterative model estimation   #######
  ####################################################
  
  library(xgboost)
  library(lme4)
  it=1
  converged=FALSE
  while(!converged && it<itmax) {
    
    #xgb
    target=rep(0,N) #target=eta-Z%*%b
    for (i in 1:N) {
      b.temp=as.matrix(bi[group[i]], nrow=q, ncol=1)
      z.temp=as.matrix(Zi.int[i,], nrow=1, ncol=q)
      target[i]= eta.est[i] - z.temp%*%b.temp
    }
    xgb.data=cov[xnam] |> data.matrix()
    xgb=xgboost(
      data=xgb.data, label=target,
      params=list(
        subsample=subsample, eta=eta, max_depth=max_depth,
        lambda=lambda, alpha=alpha
      ), 
      nrounds=nrounds,
      verbose = 0
    ) 
    f.x_ij=predict(xgb, newdata=xgb.data)
    
    #glm con mixed effects
    #glmer.data=cbind(y,group)
    #if(z.not.null) glmer.data=cbind(glmer.data, Zi)
    #glmer.data=data.frame(glmer.data) #otherwise glmer doesn't work
    glmer.data = data.frame(y,group)
    if(z.not.null) glmer.data=data.frame(glmer.data, Zi)
    glmer.fit= glmer(glmer.formula, glmer.data, family=family, offset=f.x_ij)
    
    #keep the order of the elements of the b_i
    select=c("(Intercept)",znam) #tutti i b_i da estrarre
    glmer.bi=ranef(glmer.fit)$group[select]
    
    #convergenza dei b_i
    bi.old=bi
    bi=data.frame(t(glmer.bi))
    names(bi)=lev
    diff.t=abs(bi.old-bi)
    n. diff = max (diff.t) #use the infinity (max) rule
    ind=which(diff.t==n.diff, arr.ind=T)
    n.old=abs(bi.old[ind])
    converged= n.diff/n.old <toll
    it=it+1
    all.bi[[it]]=bi
  }
  
  ###############################################
  ####	STEP 5: Output preparation       ############
  ###############################################
  
  #If non-convergence, it gives an error message
  if(!converged) {
    warning(If non-convergence, it gives an error message 'Maximum number of iterations exceeded, no convergence reached')
  }
  
  result=list(glmer.fit,xgb,bi,it,converged,all.bi,linkf,linkinv,xnam,znam,family)
  names(result)=c('glmer.model', 'xgb.model', 'rand.coef', 'n.iteration',
                  'converged','all.rand.coef','linkf','linkinv','xgb.var',
                  'random.eff.var','family')
  class(result)='gmexgb'
  result
}



summary.gmexgb=function(gm) {
  print('Mixed effects model') #summary of the mixed effects model
  print(summary(gm$glmer.model)) 
  str=ifelse(gm$converged, 'Converged', 'Did not converge')
  print(paste(str , 'after', gm$n.iteration, 'iterations')) 
}
fitted.gmexgb=function(gm, type='response', alpha=0.5) {
  #Are the fitted values from the lmer function correct?
  #I incorporated the xgb fitted values into the model
  mu=fitted(gm$glmer.model)
  #error if the response type is not among these
  allowed=c('response', 'mu', 'eta')
  msg=paste('Possible choices are ', allowed[1],', ' ,allowed[2],' and ', 
            allowed[3], sep='')
  if(sum(type==allowed)==0) stop('Type of prediction not available:',msg)
  if(type=='mu') ans=mu
  if(gm$family=='binomial' && type=='response') ans=mu>alpha
  if(type=='eta') ans=gm$linkf(mu)
  ans
}



predict.gmexgb=function(gm, newdata, group, type='response', alpha=0.5, 
                        predict.all=FALSE, re.form=NULL, newparam=NULL, 
                        terms=NULL, allow.new.levels=TRUE, na.action=na.pass,
                        random.only=FALSE) {
  
  #call the predict methods of gmexgb and xgb
  xgb.data=newdata[gm$xgb.var] |> data.matrix()
  glmer.data=data.frame(newdata[gm$random.eff.var],group)
  p1=predict(gm$glmer.model,newdata=glmer.data,newparam=newparam,re.form=re.form,
             random.only=random.only, type='link', na.action=na.action,
             allow.new.levels=allow.new.levels)			
  p2=predict(gm$xgb.model,xgb.data,predict.all=predict.all)
  
  #error if the response type is not among these
  allowed=c('response', 'mu', 'eta')
  msg=paste('Possible choices are ', allowed[1],', ' ,allowed[2],' and ', 
            allowed[3], sep='')
  if(sum(type==allowed)==0) stop('Type of prediction not available:',msg)
  
  eta=p1+p2 
  mu=gm$linkinv(eta)
  if(type=='eta') ans=eta #eta:predico eta
  if(type=='mu') ans=mu #predico la probabilit?
  if(gm$family=='binomial' && type=='response') ans=mu>alpha #Predict the answer based on alpha
  ans
}
