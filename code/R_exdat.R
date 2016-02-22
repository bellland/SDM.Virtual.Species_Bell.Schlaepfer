## based on https://pvanb.wordpress.com/2014/05/13/a-new-method-and-tool-exdet-to-evaluate-novelty-environmental-conditions/

## NT2

# Names of the reference (ref) and projection (pro) data
ref <- c("AusBio13.asc", "AusBio14.asc", "AusBio5.asc", "AusBio6.asc")
pro <- c("SaBio13.asc", "SaBio14.asc", "SaBio5.asc", "SaBio6.asc")
 
# Import the data in R using the read.asciigrid function 
# of the sp package
refdat <- prodat <- list()
for(i in 1:4){
    refdat[[i]] <- read.asciigrid(fname=ref[i])@data
    prodat[[i]] <- read.asciigrid(fname=pro[i])@data
}
refdat <- do.call(cbind, refdat)
prodat <- do.call(cbind, prodat)

# OR with raster
library(raster)
refdat <- as.matrix(stack(ref))
prodat <- as.matrix(stack(pro))

 
# Calculate the average and covariance matrix of the variables 
# in the reference set
ref.av  <- colMeans(refdat, na.rm=TRUE)
ref.cov <- var(refdat, na.rm=TRUE)
 
# Calculate the mahalanobis distance of each raster 
# cell to the environmental center of the reference 
# set for both the reference and the projection data 
# set and calculate the ratio between the two.
mah.ref    <- mahalanobis(x=refdat, center=ref.av, cov=ref.cov)
mah.pro   <- mahalanobis(x=prodat, center=ref.av, cov=ref.cov)
mah.max <- max(mah.ref[is.finite(mah.ref)])
nt2 <- as.data.frame(mah.pro / mah.max)
 
# Create and plot the raster layer
NT2 <- read.asciigrid(fname=pro[1])
NT2@data <- nt2
library(raster)
NT2rast <- raster(NT2)
plot(NT2rast, col=rainbow(100), ylim=c(-35,-20), xlim=c(15,35))




## comment by Matthew Bayly
#——————————————————————–#
# NT1 – UNIVARIATE EXTRAPOLATION
#——————————————————————–#
# UD(ij) = min(
MaxMin <- matrix(NA, nrow=2, ncol=ncol(refdat))
colnames(MaxMin) <- colnames(refdat); rownames(MaxMin) = c("max", "min")
for(i in 1:ncol(refdat)){
MaxMin[1,i] <- max(refdat[,i], na.rm=TRUE) # max value of each variable in the reference matrix
MaxMin[2,i] <- min(refdat[,i], na.rm=TRUE) # min value of each variable in the reference matrix
}
UDs <- matrix(NA, nrow(prodat), ncol(prodat))
for(j in 1:ncol(prodat)){
for(i in 1:nrow(prodat)){
UDs[i,j] <- (min(c(
(prodat[i,j] – MaxMin[2,j]),(MaxMin[1,j] – prodat[i,j]), 0), na.rm=TRUE))/
(MaxMin[1,j] – MaxMin[2,j])
}
}
UDs <- data.frame(UDs)
NT1 <- rowSums(UDs, na.rm=TRUE)

# Create and plot the raster layer NT1
mah.proR <- read.asciigrid(fname=pro[1]) # coordinates from a single grid
NT1 <- data.frame(NT1)
mah.proR@data <- NT1 # assing data values from above
library(raster)
NT1rast <- raster(mah.proR)
plot(NT1rast, col=rainbow(100), ylim=c(-35,-20), xlim=c(15,35))
