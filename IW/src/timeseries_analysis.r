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

merge_time_frames <- function(path,timesamples){
    meta  <- read_meta(path)
    files <- list.files(path=path,pattern="^run")
    df   <- data.frame()
    for (f in files[timesamples] ){
        fname =  paste(path,f,sep='/')
        frame <- read_feather(fname)
        df    <- rbind(df,frame)
    }
    return ( df )

}

###################################################
#              PLOTS + MISC                       # 
#                                                 # 
###################################################


takespec <- function(ds,s=1:30,scale=1,xlab="Freq",
                     ylab="Spec Density",title="Spectrum"){
    spec <- spectrum(ds$disp,plot=FALSE)
    plot <- plot(spec$freq[s]*scale,spec$spec[s],
                 xlab=xlab,
                 ylab=ylab)
    title(title)
}

