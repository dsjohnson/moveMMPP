
library(Matrix)
library(magrittr)
library(moveMMPP)
v <- c(85,10,5) %>% {t(./sum(.))}
M <- matrix(c(1,0.5,0,1,1,0.5,0,0.5,0.5),3,3); diag(M) <- -diag(M)
L <- diag(c(0.5, 0.8, 1)) 
MmL <- M-L

moveMMPP:::v_exp_QmLd(v, MmL)
expQ2::v_exp_Q(v,MmL, prec=1.0e-8, renorm=FALSE, check=FALSE)
v%*%expm::expm(M-L)
expm::expAtv(t(M-L), v)$eAtv

# Transpose for smoothing...
moveMMPP:::v_exp_QmLd(v, t(M))
expm::expm(M)%*%t(v)
expm::expAtv(M, v)
