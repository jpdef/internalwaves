#!/usr/bin/Rscript

# Source Libraries
library('rjson')

source('src/timeseries_analysis.r')
source('src/inversion.r')
source('src/modes.r')

#cpath   <- readline(prompt="Path : ")
args = commandArgs(trailingOnly=TRUE)

if (length(args) == 1) {
    print("Running inversion script...")
}else{
    print("Usage need config file as input")
    exit(1)
}

#Setup datapath and read metafile
inputs  <- fromJSON(file=args[1])
meta    <- read_meta(inputs$path)


#Ocean profile
stratdf <- read_feather('../config/strat.fthr')
strat   <- stratdf$strat
depth   <- seq(0,meta$depth_max,length.out=meta$depth_len)


#Sampling and Noise Params
recordstop <- 578-1*96
timesamples <- seq(1,recordstop,4)
epsi        <- 1e-3

#Read in the dataV
ds <- read_data_dir(inputs$path,timesamples,epsi)


#Make a initial parameter space
modes <- inputs$modes + 1 #zero offset in python
#omega <- inputs$freqs/3600
omega  <- c(0.0700,0.0805,0.0900)/3600
k <- generate_wavenumbers(depth,strat,omega,modes)


#Make a set of kx , ky
headings  <- seq(0,360,10)
kx <- ky <- c()
i <- 1
for (ks in k) {
        kx <- c(kx,ks*cos( pi * headings/180))
            ky <- c(ky,ks*sin( pi * headings/180))
                i <- i + 1
}

ps <- data.frame(k=rep(k,each=length(headings)),heading=headings, kx=kx,ky=ky,
                 omega=rep(omega,each=length(modes)*length(headings)),
                                  modes=rep(modes,each=length(headings)) ) 

ps <- ps[order(ps$omega,ps$modes),]



#Invert Matrix, evaluate parameters
H <- make_inv_matrix(ds,ps)

print("Inverting matrice may take time...")
system.time(approx <- tappered_least_square(H,ds$disp))

xl <- length(approx$x_hat)
ps$a  <- approx$x_hat[1:(xl/2),]
ps$b  <- approx$x_hat[(xl/2 + 1) :xl,]

library(gridExtra)
mag <- function(a,b){
        return ((a**2 + b**2))
}


par(mfrow=c(2,1))
#Plot Square Amplitudes
ggplot(ps,aes(x=kx,y=ky,size=mag(a,b),color=mag(a,b))) + 
       geom_point() + 
       scale_color_gradient(low = "blue", high = "red")


#Read in the next day
timesamples <- seq(recordstop,578,4)
ds_fut      <- read_data_dir(inputs$path,timesamples,epsi)
ds_fut      <- rbind(ds,ds_fut)

ds_fut     <- make_forecast(ds_fut,ps) 
rms_error  <- forecast_error(ds_fut)

plot(rms_error,type='l',xlab="Time Hours")


