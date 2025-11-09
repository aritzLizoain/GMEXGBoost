# Hyper parameters
eta=0.3
nrounds=100
max_depth=6
subsample=1
lambda=1
alpha_xgb=0
cbind(eta,nrounds,max_depth,subsample,lambda,alpha_xgb)

# Hyper parameters CV
load("hyper_cv.Rdata")
hyper
eta=hyper$eta
nrounds=hyper$nrounds
max_depth=hyper$max_depth
subsample=hyper$subsample
lambda=hyper$lambda
alpha_xgb=hyper$alpha_xgb
cbind(eta,nrounds,max_depth,subsample,lambda,alpha_xgb)

eta=0.1
nrounds=50
max_depth=3
subsample=0.5
lambda=1
alpha_xgb=0
cbind(eta,nrounds,max_depth,subsample,lambda,alpha_xgb)

# Simulation parameters

Model_setting=rep(c("Int","Int+Slope"),each=4)
F_eff_var=rep(c("Small","High"),times=4)
Alpha=rep(c(0.4,0.7),times=4)
Beta=rep(c(0.25,0.6),times=4)
R_eff_var=rep(c("Small","High"),each=2,times=2)
Gama2=rep(c(0.5,2,0.3,1.4),each=2)
Delta2=rep(c(0,0.5,1.4),c(4,2,2))
table2=data.frame(Model_setting,F_eff_var,Alpha,Beta,R_eff_var,Gama2,Delta2)
table2

tbl2=as_tibble(table2) |> 
  mutate(Out=vector("list",8))
tbl2

N=10
n_train=40
n_test=50

nsim=100
nsim=10

tic("Simulation")
for (i in 1:8) {
  set.seed(123)
  alpha=Alpha[i]
  beta=Beta[i]
  gama2=Gama2[i]
  delta2=Delta2[i]
  
  sim1=t(replicate(
    nsim,
    sim_pref_total(
      N,n_train,n_test,alpha,beta,gama2,delta2,
      eta=eta,nrounds=nrounds,
      max_depth=max_depth,subsample=subsample,
      lambda=lambda,alpha_xgb=alpha_xgb
    )
  ))
  # sim1
  
  tbl2$Out[[i]]=table_summary(sim1)
}
toc()

tbl2
tbl3=tbl2 |> 
  unnest(Out)
tbl3

excel_file=paste0("gmexgb_",Sys.Date(),".xlsx")
excel_file
write_xlsx(tbl3,excel_file)

tbl4 = tbl3 |> 
  group_by(Model_setting,R_eff_var,F_eff_var) |> 
  arrange(Mean_PMCR, .by_group = TRUE)
tbl4
