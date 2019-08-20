#TS_ANALYSIS
#Desc : Library for time series analysis in R
#Auth : J. DeFilippis
#Date : 7-18-2019

library(feather)
library(ggplot2)


###################################################
#               SERIES EXTRACTION                 # 
#                                                 # 
###################################################


read_meta <- function(path){
    fname = paste(path,"meta.fthr" , sep='/')
    meta <- read_feather(fname)
    return (meta)
}

depth_series <- function(path,time,xpnt){
    #Get dimensions from meta file
    meta <- read_meta(path)
    dl   <- meta$depth_len
    
    fname = paste('run-',time,'.fthr',sep="")
    fpath = paste(path,fname,sep='/')
    df    = read_feather(fpath)
    return ( df[ seq(xpnt,nrow(df),dl), ] )
}


range_series <- function(path,time,zpnt){
    #Get dimensions from meta file
    meta <- read_meta(path)
    rl   <- meta$range_len
    
    fname = paste('run-',time,'.fthr',sep="")
    fpath = paste(path,fname,sep='/')
    df    = read_feather(fpath)
    
    return (df[seq(1+ rl*(zpnt-1),rl*zpnt,1), ])
}


time_series <- function(path,zpnt,xpnt){
    #Get dimensions from meta file
    meta <- read_meta(path)
    rl   <- meta$range_len
    
    # Grab all time frames for internal wave field
    values = c()
    for (file in list.files(path=path,pattern="^run")){
        fname = paste(path,file , sep='/')
        frame <- read_feather(fname)
        slicedepth = frame[seq(1+ rl*(zpnt-1),rl*zpnt,1), ]
        slicerange = slicedepth[xpnt,]
        values <- c(values, slicerange$disp)
        
    }
    time <- seq(0,meta$time_max,length.out=meta$time_len)
    ts = data.frame(time=time,disp=values)
   
    return(ts)
}

###################################################
#          MULTI POINT EXTRACTION                 # 
#                                                 # 
###################################################


coord_to_row <- function(coord,meta){
    return ( coord[1] + meta$range_len*(coord[3]-1)) 
}

generate_random_coords <- function(num_of_points,meta){
    x <-sample(seq(1,meta$range_len),num_of_points)
    z <-sample(seq(1,meta$depth_len),num_of_points)
    coords<-list()
    for ( i in seq(1,num_of_points)){
        coords[[i]]<- c(x[i],0,z[i])
    }
    return (coords)
}

generate_coords <- function(x,z,meta){
    coords<-list()
    for ( i in seq(1,length(z))){
        for (j in seq(1,length(x))){
            coords[[i+length(z)*(j-1)]] <- c(x[j],0,z[i])
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

###################################################
#              PLOTS + MISC                       # 
#                                                 # 
###################################################


ts_compare_plot<-function(path,zpnt,xpnt1,xpnt2){
    ts1 = time_series(path,zpnt,xpnt1)
    ts2 = time_series(path,zpnt,xpnt2)
    mm =  max_min(ts1$disp,ts2$disp)
    plot(ts1,col='blue',type='l',
         xlab="Time (hours)", ylab="Displacement (m)",
         ylim= mm + mm/2) 
    lines(ts2,col='red')
    
    legend("topleft",
           legend=c(paste("X =",xpnt1/10,"km"), paste("X =",xpnt2/10,"km")),
           col=c("red", "blue"),
           lty=c(1,1),
           cex=0.7) 
    
    title(path)
}


max_min <- function(ts1,ts2){
    max <- if ( max(ts1) > max(ts2) ) max(ts1) else max(ts2)
    min <- if ( min(ts1) < min(ts2) ) min(ts1) else min(ts2)
    return (c(min,max))
}

takespec <- function(ds,s=1:30,scale=1,xlab="Freq",
                     ylab="Spec Density",title="Spectrum"){
    spec <- spectrum(ds$disp,plot=FALSE)
    plot <- plot(spec$freq[s]*scale,spec$spec[s],
                 xlab=xlab,
                 ylab=ylab)
    title(title)
}

crosscor<-function(path,zpnt,xpnt1,xpnt2){
    ts1 <- time_series(path,zpnt,xpnt1)
    ts2 <- time_series(path,zpnt,xpnt2)

    
    rho = ccf(ts1$disp,ts2$disp,plot=FALSE)
    
    plot(rho$lag,rho$acf,col='black',type='h',
         xlab="Lag", ylab="Correlation") 
    title(path)
}
