###################################################
#                  MODE FUNCTIONS                 # 
#                                                 # 
###################################################
library(pracma)
library(geigen)

cannonical_sigma <-function(depth){
    d <- max(depth)
    sigma <- 22 + 2.5*tanh(2*pi*(depth - .15*d)/d)
    return(sigma)
}


cannonical_bfrq<-function(depth){
    deriv <- fderiv(cannonical_sigma,depth)
    N <- 2*sqrt(deriv)
    return(N/3600) # return in cps
}


finite_diff_matrix <- function(N,delta){
    M <- diag(-2,N)
    M[ abs(row(M) - col(M) ) == 1 ] <- 1
    M <- M / (delta)**2
    return (M)
}


generate_vert_modes <- function(depth,strat,freq){
    r <- range(depth)
    l <- length(depth)
    delta <- (r[[2]] - r[[1]])/l
    fsq <- freq**2
    
    D2 <- finite_diff_matrix(l,delta)
    F2 <- diag(fsq,l)
    N2 <- diag(strat)
    
    M1 <- F2-N2
    return( geigen (M1,D2,symmetric = FALSE) )
}


# Note factor of 7 needs to be fixed disagreement b/t
# python modes and R modes
evaluate_mode <- function(ev,depth,mode,z){
    f <- approxfun(depth,ev$vectors[,mode])
    col <- f(z)/7
    return(col)
}


generate_wavenumbers <- function(depth,strat,omega,modes){
    wn <- c()
    for (o in unique(omega)){
        e <- generate_vert_modes(depth,strat,o)
        wnn <- o/sqrt( e$values[modes] )
        wn <- c(wn,wnn) 
   }
    return(wn)
}


generate_mode_matrix <- function(depth,strat,ps,ds){
    freqs <- unique(ps$omega)
    modes <- unique(ps$modes)
    M <- matrix(nrow=nrow(ds))
    for (f in freqs){
       
       #Eigen value/vector matrix
       ev <-  generate_vert_modes(depth,strat,f)
       
       for (m in modes) {
           mv <-  evaluate_mode(ev,depth,m,ds$z)
           #Number of wavenumbers per mode,frequency
           nk <-  nrow(ps[ps$omega == f & ps$modes == m,])
           Mm <-  matrix(mv,length(mv),nk)
           M  <-  if (all(is.na(M))) Mm else cbind(M,Mm)
       }
     }
     return (M)
}


