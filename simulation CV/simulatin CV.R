alpha=0.4
beta=0.25
gama2=0.5
delta2=0

N=10
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

# CV

Nsamp <- nrow(trainset)
Nsamp
set.seed(123)
fold_number <- sample(1:5, Nsamp, replace=TRUE)
table(fold_number)
library(tidyr)
args(gmexgb)
p <- expand_grid(
  eta=c(.1, .5, .9),
  nrounds=c(50,100,150),
  max_depth=c(3, 6, 12),
  subsample=c(0.5,0.75,1),
  lambda=1,
  alpha_xgb=0
)
p

head(trainset)

params=as.matrix(p)
n_grid=nrow(params)
n_grid
error <- matrix(0, nrow=n_grid, ncol=5)
library(progress)
pb <- progress_bar$new(total = n_grid*5)
for(i in 1:n_grid){
  for(k in 1:5){
    fold_idx <- (1:Nsamp)[fold_number == k]
    fit1=gmexgb(
      y=trainset$y[-fold_idx],
      cov=trainset[-fold_idx,],
      group = trainset$group[-fold_idx], 
      xnam=paste0("x",1:7),
      eta=params[i, 'eta'],nrounds=params[i, 'nrounds'],
      max_depth=params[i, 'max_depth'],subsample=params[i, 'subsample'],
      lambda=params[i, 'lambda'],alpha=params[i, 'alpha_xgb'],
      family='binomial'
    )
    pred_resp <-  ifelse(predict.gmexgb(fit1,trainset[fold_idx,],group=trainset$group[fold_idx],type='response')== FALSE, "0", "1") 
    error[i, k]=mean(pred_resp != trainset$y[fold_idx])
    pb$tick()
  }
}

e=rowMeans(error)
id_min=which.min(e)
id_min
p[id_min,]
params=bind_cols(p,e=e)
library(ggplot2)
ggplot(params,aes(x=eta,y=e))+geom_point()
ggplot(params,aes(x=nrounds,y=e))+geom_point()
ggplot(params,aes(x=subsample,y=e))+geom_point()
ggplot(params,aes(x=lambda,y=e))+geom_point()
ggplot(params,aes(x=alpha_xgb,y=e))+geom_point()

hyper=bind_cols(p,e=e) |> 
  dplyr::slice(id_min) |> 
  mutate(
    alpha=alpha,
    beta=beta,
    gama2=gama2,
    delta2=delta2,
    .before=1
  )
hyper

save(hyper,file = "hyper_cv.Rdata")
load("hyper_cv.Rdata")

