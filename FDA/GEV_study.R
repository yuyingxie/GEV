# TODO: Add comment
# 
# Author: xyy
###############################################################################
library(PMA)
library(MASS)
library(RSpectra)
library(msda)
library(dsda)
source("GEV.R")

args=(commandArgs(TRUE))

##args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
if(length(args)==0){
    case.id = 1  
}else{
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}

load(paste("Data/", Type, "_p", p, "_n_", n, "_id", case.id, ".Rdata", sep = ""))

Sig_est = Ome_est = matrix(0, 2 * p, 2 * p)

Ome_est[1:p, (p + 1): (2 * p)] = t(dat$X) %*% dat$Y / n
Ome_est[ (p + 1): (2 * p), 1:p] = t(dat$Y) %*% dat$X /n
Sig_est[1:p, 1:p] = t(dat$X) %*% dat$X /n
Sig_est[(p + 1): (2 * p), (p + 1): (2 * p)] = t(dat$Y) %*% dat$Y / n

perm.out <- CCA.permute(x = dat$X, z = dat$Y, typex="standard", standardize = FALSE,
        typez = "standard", nperms = 10, penaltyxs = seq(.02, .9, len = 10))
W1 = CCA(x = dat$X, z = dat$Y,  penaltyx = perm.out$bestpenaltyx, standardize = FALSE,
        v = perm.out$v.init, penaltyz = perm.out$bestpenaltyz, K = 2)

U_est = Re(Get_U(Ome_est, 2))
res = CCACV(dat$X, dat$Y, fold = 5, lambdax =seq(0.01, 0.2, length =20), 
        lambday =seq(0.01, 0.2, length = 20), k = 2)
lambda = c(rep( res$lambdaopt[1], p), rep( res$lambdaopt[2], p))
W2 = GEV_FISTA(Sig_est, U_est, lambda = lambda, diff_thre = 1e-6, max_iter = 1000)

result = list(PMA_A = Get_errorF(dat$A, W1$u), PMA_B = Get_errorF(dat$B, W1$v), 
                    GEV_A = Get_errorF(dat$A, W2[1:p, ]), GEV_B = Get_errorF(dat$B, W2[(p + 1):(2*p), ]), 
                    resCV = res)
save(result,  file = paste("./Result/",  Type, "_p", p, "_n_", n,"_id", case.id,"_Res.Rdata", sep = "")) 
