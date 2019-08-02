#TS_ANALYSIS
#Desc : Library from time series analysis in R
#Auth : J. DeFilippis
#Date : 7-18-2019

library(feather)
library(ggplot2)


depth_series <- function(path,time,xpos){
    fname = paste('run-',time,'.fthr',sep="")
    fpath = paste(path,fname,sep='/')
    df    = read_feather(fpath)
    return ( df[ seq(xpos,nrow(df),100), ] )
}


range_series <- function(path,time,zpos){
    fname = paste('run-',time,'.fthr',sep="")
    fpath = paste(path,fname,sep='/')
    df    = read_feather(fpath)
    return (df[seq(1+ 100*(zpos-1),100*zpos,1), ])
}


time_series <- function(path,zpos,xpos){
    # Grab all time frames for internal wave field
    values = c()
    for (file in list.files(path=path,pattern="^run")){
        fname = paste(path,file , sep='/')
        frame <- read_feather(fname)
        slicedepth = frame[seq(1+ 100*(zpos-1),100*zpos,1), ]
        slicerange = slicedepth[xpos,]
        values <- c(values, slicerange$disp)
        
    }

    fname = paste(path,"meta.fthr" , sep='/')
    meta <- read_feather(fname)
    ts = data.frame(time=meta$time,disp=values)
   
    return(ts)
}


ts_compare_plot<-function(path,zpos,xpos1,xpos2){
    ts1 = time_series(path,zpos,xpos1)
    ts2 = time_series(path,zpos,xpos2)
    mm =  max_min(ts1$disp,ts2$disp)
    plot(ts1,col='blue',type='l',
         xlab="Time (hours)", ylab="Displacement (m)",
         ylim= mm + mm/2) 
    lines(ts2,col='red')
    
    legend("topleft",
           legend=c(paste("X =",xpos1/10,"km"), paste("X =",xpos2/10,"km")),
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

crosscor<-function(path,zpos,xpos1,xpos2){
    ts1 <- time_series(path,zpos,xpos1)
    ts2 <- time_series(path,zpos,xpos2)

    
    rho = ccf(ts1$disp,ts2$disp,plot=FALSE)
    
    plot(rho$lag,rho$acf,col='black',type='h',
         xlab="Lag", ylab="Correlation") 
    title(path)
}
