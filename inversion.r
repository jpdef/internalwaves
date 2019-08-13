#MODEL INVERSION
#Desc : Library for model inversion in R
#Auth : J. DeFilippis
#Date : 8-08-2019


#Creates a model manifold from a set of basis functions
make_manifold <- function(samples,params,basis){
   Hv<- c()
   for (s in samples){
       row <- c()
       for (b in basis){
           row <- c(row,b(params*s))
       }
       Hv <- c(Hv,row)
   }
   H <- matrix(Hv,nrow=length(samples),byrow=TRUE)
   return (H)
}

#Makes a tridiagonal matrix with 1 on the diag and -1 
# on the elements on either side of the diagonal
make_smoothness_matrix <- function(n){
    m <- diag(1,n)
    m[ abs( row(m) - col(m)) == 1 ] <- -1
    return ( m )
}

#Makes a weight matrix based on the variances 
# of a random variable that has a guassian distrubution
make_weight_matrix <- function(n){
    Mv <- c()
    r <- rnorm(n,mean=1)
    for (i in seq(1,n) ){
        Mv <- c(Mv,r[i]*r)
    }
    M <- matrix(Mv,nrow=n,byrow=TRUE)
    
    return (M)
}

#Tappered least square method wunsch
tappered_least_square_w <- function(samps,obs,params,basis_fn){
    H <- make_manifold(samps,params,basis_fn)
    W <- diag(length(obs))
    S <- make_smoothness_matrix(length(basis_fn)*length(params))
    M1 <- H %*% S %*% t(H)
    x_approx <- S %*% t(H) %*% solve( M1 + W ) %*% obs
    y_approx <- H%*%x_approx
    n_approx <- obs-y_approx
    return( list( x=x_approx,y=y_approx,n=n_approx ) )
    
}

#Tappered least square method cornuelle
tappered_least_square_c <- function(samps,obs,params,basis_fn){
    H <- make_manifold(samps,params,basis_fn)
    R <- diag(length(obs))
    Q <- diag(length(basis_fn)*length(params))
    R_inv <- solve(R)
    Q_inv <- solve(Q)
    x_approx <- solve(t(H) %*% R_inv %*% H + Q_inv) %*% t(H)%*% R_inv %*% obs
    y_approx <- H%*%x_approx
    n_approx <- obs-y_approx
    return( list( y=obs,x_hat=x_approx,y_hat=y_approx,n_hat=n_approx ) )
    
}

#Tappered least square method cornuelle
tappered_least_square <- function(H,obs){
    R <- 0.1*diag(nrow(H))
    Q <- 0.1*diag(ncol(H))
    R_inv <- solve(R)
    Q_inv <- solve(Q)
    x_approx <- solve(t(H) %*% R_inv %*% H + Q_inv) %*% t(H)%*% R_inv %*% obs
    y_approx <- H%*%x_approx
    n_approx <- obs-y_approx
    return( list( y=obs,x_hat=x_approx,y_hat=y_approx,n_hat=n_approx ) )
    
}

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
    #return( geigen (M1,D2,symmetric = FALSE) )
    return( eigen (solve(D2) %*% M1) )
}

evaluate_modes <- function(ev,depth,modes,z){
    fm <- c()
    for (m in modes){
        f <- approxfun(depth,ev$vectors[,m])
        fm <- c(fm,f)
    }
    vm <- c()
    for (f in fm){
        col <- f(z)
        vm <- c(vm,col)
    }
    return(matrix(vm,nrow=length(z)))
}

generate_wavenumbers <- function(depth,strat,ps){
    modes <- unique(ps$modes)
    wn <- c()
    for (o in unique(ps$omega)){
        e <- generate_vert_modes(depth,strat,o)
        wnn <- o/sqrt( e$values[modes] )
        wn <- c(wn,wnn) 
    }
    return(wn)
}


generate_mode_matrix <- function(depth,strat,ps,ds){
    freqs <- unique(ps$omega)
    modes <- unique(ps$modes)
    for (f in freqs){
       ev <-  generate_vert_modes(depth,strat,f)
       Mn <-  evaluate_modes(ev,depth,modes,ds$z)
       M  <- if( exists('M') ) cbind(M,Mn) else Mn
     }
     return (M)
}

###################################################
#                  PLOT FUNCTIONS                 # 
#                                                 # 
###################################################

plot_ls_fit <- function(fit,xaxis,xlab,ylab){
    plot(xaxis,fit$y,col='black',type='l',
              xlab=xlab,
              ylab=ylab)
    lines(xaxis,fit$y_hat,col='red')
    title("Plot of Fit")
}

plot_ls_res <- function(fit,xaxis,xlab,ylab){
    plot(xaxis,fit$n_hat,col='black',type='p',
              xlab=xlab,
              ylab=ylab)
    title("Plot of Residuals")
}

hist_ls_res <- function(fit){
    n = fit$n_hat
    hist(n,xlab="Residuals",main="Histogram of Residuals")
}

mag <- function(a,b){
    return (sqrt(a**2 + b**2))
}

plot_ls_par <- function(fit,params,parlab,par_true,xlab,ylab){
    colors <- list('blue','red','green')
    for (m in unique(ps$modes) ){
        s = ps$mode == m
        if (m == 1){
            y <- ps[s,parlab]
            plot(3600*ps$omega[s],y,
                 #ylim=10*c(min(y),max(y)),
                 ylim=2*c(min(fit$y),max(fit$y)),
                 xlab=xlab,
                 ylab=ylab,
                 col=colors[[m]])
        }
        points(3600*ps$omega[s],ps[s,parlab],col=colors[[m]])

    }

    points(3600*ps$omega[s],par_true,pch=3)
    legend("topright",legend=c("M1","M2","M3","Model"),
            col=c('blue','red','green','black'),
            lty=c(rep(1,3),3),
            cex=0.6)
    title(paste("Plot of",parlab))

}


plot_ls_mag <- function(fit,par_true,xlab,ylab){
    colors <- list('blue','red','green')
    for (m in unique(ps$modes) ){
        s = ps$mode == m
        if (m == 1){
            y <- mag(ps$a[s],ps$b[s])
            plot(3600*ps$omega[s],y,
                 ylim=c(0,2*max(y)),
                 xlab=xlab,
                 ylab=ylab,
                 col=colors[[m]])
        }
        points(3600*ps$omega[s],mag(ps$a[s],ps$b[s]),col=colors[[m]])

    }

    points(3600*ps$omega[s],par_true,pch=3)
    legend("topright",legend=c("M1","M2","M3","Model"),
            col=c('blue','red','green','black'),
            lty=c(rep(1,3),3),
            cex=0.9)
    title("Plot of Parameters")
}

