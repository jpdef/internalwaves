#TS_ANALYSIS
#Desc : Library for time series analysis in R
#Auth : J. DeFilippis
#Date : 7-18-2019

library(feather)
library(ggplot2)


###################################################
#          MULTI POINT EXTRACTION                 # 
#                                                 # 
###################################################


read_meta <- function(path){
    fname = paste(path,"meta.fthr" , sep='/')
    meta <- read_feather(fname)
    return (meta)
}

coord_to_row <- function(coord,meta){
    l <- meta$range_len
    return ( coord[1] + l*(coord[2]-1) + l*l*(coord[3]-1)) 
}

generate_random_coords <- function(num_of_points,meta){
    x <-sample(seq(1,meta$range_len),num_of_points)
    y <-sample(seq(1,meta$range_len),num_of_points)
    z <-sample(seq(1,meta$depth_len),num_of_points)
    coords<-list()
    for ( i in seq(1,num_of_points)){
        coords[[i]]<- c(x[i],y[i],z[i])
    }
    return (coords)
}

generate_coords <- function(x,y,z,meta){
    coords<-list()
    for ( k in seq(1,length(z))){
        for (j in seq(1,length(y))){
            for (i in seq(1,length(x))){
                I <- length(x)
                J <- length(y)
                index <- i + J*(j-1) + I*J*(k-1) 
                coords[[index]] <- c(x[i],y[j],z[k])
            }
        }
    }
    return (coords)
}


multidim_samples <- function(path,coords,timesamples){
    meta  <- read_meta(path)
    files <- list.files(path=path,pattern="^run")
    t <- 1
    df   <- data.frame()
    for (f in files[timesamples] ){
        fname =  paste(path,f,sep='/')
        frame <- read_feather(fname)
        tt  <- timesamples[t]*(meta$time_max)/(meta$time_len-1)
        frame <- cbind(frame,time=tt)
        for (coord in coords ) {
            row   <- coord_to_row(coord,meta)
            df    <- rbind(df,frame[row,])
        }
        t <- t + 1
    }
    return ( df )

}


subsample <- function(ds,xsamples=c(),ysamples=c(),zsamples=c()){
        if (length(xsamples) != 0 ){
           xf <- unique(ds$x)[xsamples]
           ds <- ds[ds$x %in% xf,]
        }
        
        if (length(ysamples) != 0 ){
           yf <- unique(ds$y)[ysamples]
           ds <- ds[ds$y %in% yf,]
            
        }
        
        if (length(zsamples) != 0 ){
           zf <- unique(ds$z)[zsamples]
           ds <- ds[ds$z %in% zf,]
            
        }
    return(ds)
}

dipsample <- function(ds,xsamples=c(),ysamples=c(),zsamples=c(),center=c()){
        if (length(xsamples) != 0 ){
           xf <- unique(ds$x)[c(xsamples,center[1])]
           yc <- unique(ds$y)[center[2]]
           ds1 <- ds[ds$x %in% xf & ds$y == yc,]
        }
        
        if (length(ysamples) != 0 ){
           yf <- unique(ds$y)[c(ysamples,center[2])]
           xc <- unique(ds$y)[center[1]]
           ds2 <- ds[ds$y %in% yf & ds$x == xc,]
            
        }
        
        ds <- rbind(ds1,ds2)
    
        if (length(zsamples) != 0 ){
           zf <- unique(ds$z)[zsamples]
           ds <- ds[ds$z %in% zf,]
            
        }
    return(ds)
}




merge_time_frames <-function(path,timesamples=c(),xsamples=c(),ysamples=c(),zsamples=c(),center=c()){
    meta  <- read_meta(path)
    files <- list.files(path=path,pattern="^run")
    df   <- data.frame()
    for (f in files[timesamples] ){
        fname =  paste(path,f,sep='/')
        frame <- read_feather(fname)
        if (length(center)==0){
            frame <- subsample(frame,xsamples,ysamples,zsamples)
        }else{
            frame <-dipsample(frame,xsamples,ysamples,zsamples,center)
            
        }
        df    <- rbind(df,frame)
    }
    return ( df )

}


read_data_dir <- function(path,epsi,timesamples=c(),xsamples=c(),ysamples=c(),zsamples=c(),center=c()){
    #Read in the data
    ds <- merge_time_frames(path,timesamples,xsamples,ysamples,zsamples,center)
    
    #Select only on correct variables
    ds <- subset(ds, select=c(x,y,z,t,d))
    
    #Add some noise
    noise <- epsi*rnorm(length(ds$d)) 
    ds$d <- ds$d + noise
    return (ds)
}


###################################################
#              PLOTS + MISC                       # 
#                                                 # 
###################################################

plot_ts <- function(ds,xi,yi,zi){
    xs <- unique(ds$x)[xi]
    ys <- unique(ds$y)[yi]
    zs <- unique(ds$z)[zi]
    ts <- ds[ ds$x==xs & ds$y == ys & ds$z == zs ,]
    plot(ts$t,ts$d,type='l')
}

takespec <- function(ds,s=1:30,scale=1,xlab="Freq",
                     ylab="Spec Density",title="Spectrum"){
    spec <- spectrum(ds$disp,plot=FALSE)
    plot <- plot(spec$freq[s]*scale,spec$spec[s],
                 xlab=xlab,
                 ylab=ylab)
    title(title)
}


#Plot a time series at a specific spacial points
plot_ts <- function(ds,xi,yi,zi){
    xs <- unique(ds$x)[xi]
    ys <- unique(ds$y)[yi]
    zs <- unique(ds$z)[zi]
    ts <- ds[ ds$x==xs & ds$y == ys & ds$z == zs ,]
    plot(ts$t,ts$d,type='l')
}