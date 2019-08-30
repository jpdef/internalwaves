#MODEL INVERSION
#Desc : Library for model inversion in R
#Auth : J. DeFilippis
#Date : 8-08-2019

#Source libraries
source('modes.r')


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


#Creates inverse matrix codifing all the basis functions
make_inv_matrix <- function(ds,ps){
    
    #Sinsoid Basis
    real <- cos(2*pi*(ds$x%*%t(ps$kx) + ds$y%*%t(ps$ky) - ds$t%*%t(ps$omega)) )
    img  <- sin(2*pi*(ds$x%*%t(ps$kx) + ds$y%*%t(ps$ky) - ds$t%*%t(ps$omega)) )
    SB <- cbind(real,img)
    
    #Modal Basis
    V <- generate_mode_matrix(depth,strat,ps,ds)
    MB <- cbind(V,V) 
    
    return (MB*SB)
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
tappered_least_square <- function(H,obs){
    R_inv <- 1*diag(nrow(H))
    Q_inv <- 1*diag(ncol(H))
    GM    <- t(H) %*% R_inv %*% H + Q_inv
    x_approx <- solve(GM) %*% t(H)%*% R_inv %*% obs
    y_approx <- H%*%x_approx
    n_approx <- obs-y_approx
    return( list( y=obs,x_hat=x_approx,y_hat=y_approx,n_hat=n_approx ) )
    
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
    for (m in unique(params$modes) ){
        s = params$mode == m
        A <- sum( abs(params[s,parlab]))
        if (m == 1){
            y <- params[s,parlab]
            plot(3600*params$omega[s],y/A,
                 #ylim=10*c(min(y),max(y)),
                 #ylim=2*c(min(fit$y),max(fit$y)),
                 ylim=c(-0.35,0.35),
                 xlab=xlab,
                 ylab=ylab,
                 col=colors[[m]])
        }
        points(3600*params$omega[s],params[s,parlab]/A,col=colors[[m]])

    }
    A <- sum(abs(par_true))
    points(3600*params$omega[s],par_true/A,pch=3)
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
            cex=0.6)
    title("Plot of Parameters")
}

