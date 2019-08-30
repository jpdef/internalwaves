# Source Libraries
library('rjson')

source('src/timeseries_analysis.r')
source('src/inversion.r')

print("Running inversion script...")

#Setup datapath and read metafile
path    <- '../data/simple/'
meta    <- read_meta(path)
inputs  <- fromJSON(file='config/run1.json')

#Ocean profile
stratdf <- read_feather('config/strat.fthr')
strat   <- stratdf$strat
depth   <- seq(0,meta$depth_max,length.out=meta$depth_len)

#Read in the data
timesamples <- seq(0,2000,2)
epsi <- 1e-3
ds <- merge_time_frames(path,timesamples)

#Add some noise
noise <- epsi*rnorm(length(ds$disp)) 
ds$disp <- ds$disp + noise
colnames(ds) <- c("x","z","disp","time")


#Make a initial parameter space
modes <- inputs$modes + 1 #zero offset in python
omega <- inputs$freqs/3600
k <- generate_wavenumbers(depth,strat,omega,modes)
ps <- data.frame(a=rep(0,length(omega)),omega=omega,k = k, modes=modes)


#Invert Matrix, evaluate parameters
H <- make_inv_matrix(ds,ps,depth,strat)
approx <- tappered_least_square(H,ds$disp)
xl <- length(approx$x_hat)
ps$a  <- approx$x_hat[1:(xl/2),]
ps$b  <- approx$x_hat[(xl/2 + 1) :xl,]

#Print out parameter space to user, compute error
amps <-  inputs$amps / length(inputs$modes)
amps <-  rep(amps,length(modes)) 
ps$atrue <- amps*cos(pi*inputs$phase/180)
ps$btrue <- amps*sin(pi*inputs$phase/180)

print(ps)
error <-  abs ( sum (ps$a - ps$atrue) ) + abs (sum (ps$b + ps$btrue) ) 
error <- error/(length(ps$atrue)+length(ps$btrue))
print(paste ( "Model Avg Error : " , format(error,digits=2) ) )
print(paste ( "Number of Samples : " , format(nrow(ds),digits=2) ) )
