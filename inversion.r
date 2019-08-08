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
tappered_least_square <- function(samps,obs,params,basis_fn){
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
    hist(n,breaks=length(n)/10,xlab="Residuals",main="Histogram of Residuals")
}

plot_ls_par <- function(fit,xaxes,xlab,ylab,xtrue){
    colors <- list('green','blue','red')
    xlim <-c(min(xaxes[[1]]),max(xaxes[[1]]))
    ylim <-c(min(fit$x_hat),max(fit$x_hat))
    if (missing(xtrue)) {
        plot(0,xlim=xlim,ylim=ylim,type='h',xlab=xlab,ylab=ylab)
    }else{
        plot(xtrue[[1]],xtrue[[2]],xlim=xlim,ylim=ylim,type='h',xlab=xlab,ylab=ylab)
    }
    
    p <- 0
    i <- 1
    for (xaxis in xaxes){
        l <- length(xaxis) + p
        points(xaxis,fit$x_hat[seq(p+1,l)],col=colors[[i]])
        p <- l
        i <- i+1
    }
    title("Plot of Parameters")
}

