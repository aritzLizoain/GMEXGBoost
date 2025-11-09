# Hyperparameters
eta=0.3
nrounds=100
max_depth=6
subsample=1
lambda=1
alpha_xgb=0

# parameters
N=10
n=40
alpha=0.4
beta=0.25
gama2=0.5
delta2=0

set.seed(123)

x1=runif(n,-1,1)
x2=runif(n,-1,1)
x3=rweibull(n,3)
x4=runif(n,-3,3)
x5=runif(n,-6,6)
x6=runif(n,-5,5)
x7=runif(n,-4,4)

tr=tree(x4,x5,x6)
fix=fixed(x1,x2,x3,x4,x5,x6,x7,alpha,beta)
library(rpart)
model=rpart(fix~x1+x2+x3+x4+x5+x6+x7)
model
par(xpd = NA) # otherwise on some devices the text is clipped
plot(model)
text(model)

b0=rnorm(1,0,sd=sqrt(gama2))
b1=rnorm(1,0,sd=sqrt(delta2))
rand=b0+b1*x1

eta=fix+rand
mu=plogis(eta)
y=rbinom(n,size=1,prob=mu)

n_train=40
n_test=50

set.seed(123)
s1=sim_gmexgb(N,n_train,n_test,alpha,beta,gama2,delta2)
head(s1)
tail(s1)

s1 |> count(group)
s1 |> count(group,train_test)
s1 |> count(train_test)
s1 |> count(group,b0,b1)
s1 |> count(group,train_test,b0,b1)
s1 |> count(group,y)

boxplot(mu~group,data=s1)

trainset=s1 |> 
  filter(train_test=="train")
testset=s1 |> 
  filter(train_test=="test")
#######################glmm xgb
library(tictoc)
tic("model fit1")
fit1=gmexgb(y=trainset$y, cov=trainset,
            group = trainset$group, 
            xnam=paste0("x",1:7),
            family='binomial') 
toc()

fit1
############################
library(caret)
pred_resp=predict.gmexgb(fit1,testset,group=testset$group,type='response')
table(pred_resp)
yhat=pred_resp |> as.integer() |> as.factor()
table(yhat)
y=testset$y
table(y)
table(yhat,y)
confusionMatrix(yhat,y)

pmcr=mean(yhat != y)
pmcr

############# roc
pred_mu=predict.gmexgb(fit1,testset,group=testset$group,type='mu')
mu=testset$mu
plot(mu,pred_mu)
pmad=mean(abs(mu-pred_mu))
pmad

library(ROCR)
ROCR_pred_test <- prediction(pred_mu,testset$y)
ROCR_perf_test <- performance(ROCR_pred_test,'tpr','fpr')
plot(ROCR_perf_test,colorize=TRUE,print.cutoffs.at=seq(0,1,by=0.1), main = "ROC Curve for Logistic Regression Model")
perf <- performance(ROCR_pred_test ,"auc")
(auc <- as.numeric(perf@y.values))

set.seed(123)
s0=sim_pref_total(N,n_train,n_test,alpha,beta,gama2,delta2)
s0
s2=sim_pref(N,n_train,n_test,alpha,beta,gama2,delta2)
s2

nsim=100
nsim=10
set.seed(123)
s3=t(replicate(
  nsim,
  sim_pref(N,n_train,n_test,alpha,beta,gama2,delta2)
))
s3

s4=na.omit(s3)
colMeans(s4)
matrixStats::colVars(s4)

nsim=100
nsim=10
set.seed(123)
s5=t(replicate(
  nsim,
  sim_pref_total(
    N,n_train,n_test,alpha,beta,gama2,delta2,
    eta=eta,nrounds=nrounds,
    max_depth=max_depth,subsample=subsample,
    lambda=lambda,alpha=alpha_xgb
  )
))
s5

m1=colMeans(s5,na.rm = T)
v1=matrixStats::colVars(s5,na.rm = T)

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

table_summary(s5)
