library(ggplot2)
library(spatstat)
library(GET)


### Randomly permuting observed times across spatial locations, producing
### a list of space-time point patterns, the observed pattern being the first

permute.times <- function(X, n.perm=999){
  # X ............ point pattern of spatial locations with time coordinates as marks
  # n.perm ....... number of random permutations of times
  
  if (!is.numeric(n.perm)){return(warning("n.perm must be a natural number."))}
  
  Xspace <- unmark(X)
  
  output <- list()
  output[[1]] <- X
  
  for (i in 1:n.perm){
    output[[i+1]] <- X %mark% sample(marks(X))
  }
  
  return(output)
}



### Compute the test statistics: S, S.space, S.time

S.stats <- function(X, sigma.t=NULL, sigma.s=NULL, times=NULL, dimyx=c(50,50), start.time=NULL, end.time=NULL){
  # X .......... marked point pattern, marks contain the temporal coordinates of observed events
  # sigma.t .... bandwidth for the temporal kernel
  # sigma.s .... bandwidth parameter for the spatial kernel
  # times ...... vector of values t_i for which to compute \hat{\lambda}(.,t_i), assumed to for a regular grid
  # dimyx ...... pixel array dimensions in the spatial domain
  # start.time, end.time ... endpoints of T, used for computing temporal edge-correction factors
  
  if (!is.numeric(times)){return(warning("times must be a numeric vector."))}
  if (!is.numeric(dimyx) | (length(dimyx)!=2)){return(warning("dimyx must be a numeric vector of length 2."))}
  
  if (is.null(start.time)){start.time <- min(X$marks)}
  if (is.null(end.time)){end.time <- max(X$marks)}
  
  # Spatial projection process
  Xspace <- unmark(X)
  
  # Output, empty so far
  est3d <- array(dim=c(dimyx[[1]],dimyx[[2]],length(times)))
  
  
  ### Non-separable estimate ###
  ##############################
  
  nonsep <- array(dim=c(dimyx[[1]],dimyx[[2]],length(times)))
  
  # Edge-correction factors, temporal part
  edgewt.t <- pnorm((end.time-X$marks)/sigma.t) - pnorm((start.time-X$marks)/sigma.t)
  # use 1/edgewt.t !
  
  for (i in 1:length(times)){
    t1 <- times[i]
    
    # Weights for the kernel estimation in 2D
    weights1 <- ker.time(si=X$marks, t=t1, sigma.t=sigma.t)/edgewt.t
    
    # Kernel estimation in 2D with weights
    dens <- density.ppp(Xspace, kernel="gaussian", weights=weights1, sigma=sigma.s, diggle=TRUE, dimyx=dimyx)
    nonsep[,,i] <- dens$v
  }
  nonsep <- pmax(nonsep,0)
  
  
  ### Separable estimate ###
  ##########################
  
  # Kernel estimation in 2D
  aux <- density.ppp(Xspace, kernel="gaussian", sigma=sigma.s, diggle=TRUE, dimyx=dimyx)
  est2d <- aux$v
  
  # Kernel estimation in 1D
  est1d <- rep(NA, times=length(times))
  for (i in 1:length(times)){
    t1 <- times[i]
    weights1 <- ker.time(si=X$marks, t=t1, sigma.t=sigma.t)/edgewt.t
    est1d[i] <- sum(weights1)
  }
  
  sep <- pmax((est2d %o% est1d)/X$n,0)
  
  
  ### Computing the test statistic values ###
  ###########################################
  
  est3d <- nonsep/sep
  est3d[is.infinite(est3d)] <- NA
  
  sum.narm <- function(x){
    return(sum(x,na.rm=TRUE))
  }
  
  S.space <- apply(est3d, 1:2, sum)*(times[2]-times[1])
  S.time <- apply(est3d, 3, sum.narm)*aux$xstep*aux$ystep
  
  
  ### Handle possible NA values due to non-rectangular observation window ###
  ###########################################################################
  
  Wpix <- as.mask(X$window, dimyx=dimyx)
  
  coordinates.2D <- get.coordinates.2D(Wpix=Wpix, dimyx=dimyx, S.space.sample.values=S.space)
  aux.vec <- as.vector(S.space)
  S.space <- aux.vec[!is.na(aux.vec)]
  
  coordinates.3D <- get.coordinates.3D(Wpix=Wpix, dimyx=dimyx, times=times, S.sample.values=est3d)
  aux.vec <- as.vector(est3d)
  S <- aux.vec[!is.na(aux.vec)]
  
  
  ### Auxiliary values for plotting ###
  #####################################
  
  width <- Wpix$xstep
  height <- Wpix$ystep
  
  
  return(list(S=S, S.space=S.space, S.time=S.time,
              coordinates.2D=coordinates.2D, coordinates.3D=coordinates.3D,
              width=width, height=height))
}



### Auxiliary kernel for estimating temporal intensity function

ker.time <- function(si, t, sigma.t){
  # Gaussian kernel
  return(dnorm(x=si, mean=t, sd=sigma.t))
}



### Auxiliary function determining pixel positions for plotting S.space

get.coordinates.2D <- function(Wpix, dimyx, S.space.sample.values){
  coordinates2D <- expand.grid(Wpix$yrow, Wpix$xcol)
  names(coordinates2D) <- c("y", "x")
  Y <- !is.na(S.space.sample.values)
  coordinates2D <- coordinates2D[(Y),]
  return(coordinates2D)
}



### Auxiliary function determining pixel positions for plotting S

get.coordinates.3D <- function(Wpix, dimyx, times, S.sample.values){
  coordinates3D <- expand.grid(Wpix$yrow, Wpix$xcol, times)
  names(coordinates3D) <- c("y", "x", "t")
  Y <- !is.na(S.sample.values[,,1])
  coordinates3D <- coordinates3D[(Y),]
  return(coordinates3D)
}



### Plotting the 3D global envelope test output,
### either the observed function, or upper or lower envelope with significance regions

plot.global_envelope.3D <- function(S.GET, coordinates.3D, what=c("obs", "lo", "hi")){
  # S.GET ... output of the function global_envelope_test for the S statistic
  # coordinates.3D ... output of the function get.coordinates.3D
  what <- match.arg(what)
  
  names(coordinates.3D) <- c("y", "x", "t")
  
  df <- cbind(as.data.frame(S.GET), coordinates.3D,
              signif.lo = S.GET$obs < S.GET$lo,
              signif.hi = S.GET$obs > S.GET$hi)
  
  switch(what,
         obs = {
           # Observed function
           p <- ggplot(df) + geom_tile(aes(x=x, y=y, fill=obs)) +
             facet_wrap(vars(t), ncol=5) + coord_fixed(ratio=1)
         },
         lo = {
           # Lower envelope with significant regions
           p <- ggplot(df) + geom_tile(aes(x=x, y=y, fill=lo)) +
             facet_wrap(vars(t), ncol=5) + coord_fixed(ratio=1) +
             geom_tile(data=df[df$signif.lo,], aes(x=x, y=y), fill="red")
         },
         hi = {
           # Upper envelope with significant regions
           p <- ggplot(df) + geom_tile(aes(x=x, y=y, fill=hi)) +
             facet_wrap(vars(t), ncol=5) + coord_fixed(ratio=1) +
             geom_tile(data=df[df$signif.hi,], aes(x=x, y=y), fill="red")
         })
  return(p)
}

### Vizualization tools for investigating the first-order separability hypothesis
### for space-time point patterns.

### Implemented by JiĹ™Ă­ DvoĹ™Ăˇk

### Version 28.3.2021

### Codes accompany the paper "Testing the first-order separability hypothesis for spatio-temporal point patterns"
### by M. Ghorbani, N. Vafaei, J. DvoĹ™Ăˇk and M. MyllymĂ¤ki,
### published in the Computational Statistics and Data Analysis journal 161, 107245.
### https://doi.org/10.1016/j.csda.2021.107245



### Vizualization based on 3D kernel estimates of intensity function ###
########################################################################

intensitySTPP <- function(X, sigma.t=NULL, sigma.s=NULL, times=NULL, start.time=NULL, end.time=NULL){
  # X .......... marked point pattern, marks contain the temporal coordinates of observed events
  # sigma.t .... bandwidth for the temporal kernel
  # sigma.s .... bandwidth parameter for the spatial kernel
  # times ...... vector of values t_i for which to compute \hat{\lambda}(.,t_i)
  # start.time, end.time ... endpoints of T, used for computing temporal edge-correction factors
  
  if (!is.numeric(times)){return(warning("times must be a numeric vector."))}
  
  if (is.null(start.time)){start.time <- min(X$marks)}
  if (is.null(end.time)){end.time <- max(X$marks)}
  
  
  ### Estimate ST intensity function by 3d kernel smoothing ###
  #############################################################
  
  est3d <- intensityST.nonsep.grid(X=X, sigma.t=sigma.t, sigma.s=sigma.s, times=times, start.time=start.time, end.time=end.time)
  
  
  ### Estimating ST intensity function in a product form ###
  ##########################################################
  
  est.sep <- intensityST.sep.grid(X=X, sigma.t=sigma.t, sigma.s=sigma.s, times=times, start.time=start.time, end.time=end.time)
  
  
  ### Standardized differences ###
  ################################
  
  est.std <- list()
  for (i in 1:length(times)){
    a <- est3d[[i]]
    b <- est.sep[[i]]
    est.std[[i]] <- eval.im((a-b)/b)
  }
  
  return(list(intensities.3d=est3d, intensities.sep=est.sep, intensities.std=est.std, times=times, sigma.s=sigma.s, sigma.t=sigma.t))
}

plot.intensitySTPP <- function(x, same.range=TRUE){
  # x           ... output of function "IntensitySTPP"
  # same.range  ... logical, should the same color scale be used for all subplots?
  
  n.patterns <- length(x$intensities.3d)
  
  if ((n.patterns==1)){
    stop("Interactive plotting not available for a single subpattern, switching to static plotting instead.")
  }
  
  # Plotting, either with the same color scale for all subpatterns
  # or with the color scale fully exploited in each subpattern
  
  ax <- list(
    visible = FALSE
  )
  ay <- list(
    visible = FALSE,
    scaleanchor = "x"
  )
  
  # ST intensity function estimated by 3d kernel smoothing
  aval1 <- list()
  for (k in 1:n.patterns){
    aval1[[k]] <- list(visible = FALSE, z=x$intensities.3d[[k]]$v, x=x$intensities.3d[[k]]$xcol, y=x$intensities.3d[[k]]$yrow)
  }
  aval1[[1]]$visible = TRUE
  
  # ST intensity function estimated in a product form
  aval2 <- list()
  for (k in 1:n.patterns){
    aval2[[k]] <- list(visible = FALSE, z=x$intensities.sep[[k]]$v, x=x$intensities.sep[[k]]$xcol, y=x$intensities.sep[[k]]$yrow)
  }
  aval2[[1]]$visible = TRUE
  
  # ST intensity function - standardized differences
  aval3 <- list()
  for (k in 1:n.patterns){
    aval3[[k]] <- list(visible = FALSE, z=x$intensities.std[[k]]$v, x=x$intensities.std[[k]]$xcol, y=x$intensities.std[[k]]$yrow)
  }
  aval3[[1]]$visible = TRUE
  
  
  if (!same.range){
    steps <- list()
    p <- p1 <- p2 <- p3 <- plot_ly()
    for (i in 1:n.patterns) {
      zmax <- max(max(aval1[[i]]$z, na.rm=TRUE),max(aval2[[i]]$z, na.rm=TRUE), na.rm=TRUE)
      zmin <- min(min(aval1[[i]]$z, na.rm=TRUE),min(aval2[[i]]$z, na.rm=TRUE), na.rm=TRUE)
      
      p1 <- add_heatmap(p1, z=aval1[[i]]$z, x=aval1[[i]]$x,  y=aval1[[i]]$y, visible = aval1[[i]]$visible, hoverinfo="none", 
                        zmax=zmax, zmin=zmin, colorbar=list(title = "Left, middle", x = 1, xanchor="left", y = 1, yanchor="top", thicknessmode="fraction", thickness=0.05, lenmode="fraction", len=0.5))
      p1 <- p1 %>% layout(xaxis = ax, yaxis = ay, title=list(text=""), margin=list(l=5,t=25,b=5))
      
      p2 <- add_heatmap(p2, z=aval2[[i]]$z, x=aval2[[i]]$x,  y=aval2[[i]]$y, visible = aval2[[i]]$visible, hoverinfo="none",
                        zmax=zmax, zmin=zmin, colorbar=list(title = "Left, middle", x = 1, xanchor="left", y = 1, yanchor="top", thicknessmode="fraction", thickness=0.05, lenmode="fraction", len=0.5))
      p2 <- p2 %>% layout(xaxis = ax, yaxis = ay, title=list(text=""), margin=list(l=5,t=25,b=5))
      
      p3 <- add_heatmap(p3, z=aval3[[i]]$z, x=aval3[[i]]$x,  y=aval3[[i]]$y, visible = aval3[[i]]$visible, hoverinfo="none", 
                        colorbar=list(title = "Right", x = 1, xanchor="left", y = 0.5, yanchor="top", thicknessmode="fraction", thickness=0.05, lenmode="fraction", len=0.5))
      p3 <- p3 %>% layout(xaxis = ax, yaxis = ay, title=list(text=""), margin=list(l=5,t=25,b=5))
      
      step <- list(args = list('visible', rep(FALSE, length(aval1))), label=x$times[i], method = 'restyle')
      step$args[[2]][i] = TRUE  
      steps[[i]] = step 
    }  
    p <- subplot(p1,p2,p3, margin=0.02) %>%
      layout(sliders = list(list(currentvalue = list(prefix = "Time coordinate: "), steps = steps)))
    return(p)
  }
  
  if (same.range){
    max.valsA <- max.valsB <- max.valsC <- rep(NA, times=n.patterns)
    min.valsA <- min.valsB <- min.valsC <- rep(NA, times=n.patterns)
    for (k in 1:n.patterns){
      max.valsA[k] <- max(aval1[[k]]$z, na.rm=TRUE)
      max.valsB[k] <- max(aval2[[k]]$z, na.rm=TRUE)
      max.valsC[k] <- max(aval3[[k]]$z, na.rm=TRUE)
      min.valsA[k] <- min(aval1[[k]]$z, na.rm=TRUE)
      min.valsB[k] <- min(aval2[[k]]$z, na.rm=TRUE)
      min.valsC[k] <- min(aval3[[k]]$z, na.rm=TRUE)
    }
    max.AB <- max(c(max.valsA,max.valsB))
    min.AB <- min(c(min.valsA,min.valsB))
    max.C <- max(max.valsC)
    min.C <- min(min.valsC)
    
    steps <- list()
    p <- p1 <- p2 <- p3 <- plot_ly()
    for (i in 1:n.patterns) {
      p1 <- add_heatmap(p1, z=aval1[[i]]$z, x=aval1[[i]]$x,  y=aval1[[i]]$y, visible = aval1[[i]]$visible, hoverinfo="none", 
                        zmax=max.AB, zmin=min.AB, colorbar=list(title = "Left, middle", x = 1, xanchor="left", y = 1, yanchor="top", thicknessmode="fraction", thickness=0.05, lenmode="fraction", len=0.5))
      p1 <- p1 %>% layout(xaxis = ax, yaxis = ay, title=list(text=""), margin=list(l=5,t=25,b=5))
      
      p2 <- add_heatmap(p2, z=aval2[[i]]$z, x=aval2[[i]]$x,  y=aval2[[i]]$y, visible = aval2[[i]]$visible, hoverinfo="none",
                        zmax=max.AB, zmin=min.AB, colorbar=list(title = "Left, middle", x = 1, xanchor="left", y = 1, yanchor="top", thicknessmode="fraction", thickness=0.05, lenmode="fraction", len=0.5))
      p2 <- p2 %>% layout(xaxis = ax, yaxis = ay, title=list(text=""), margin=list(l=5,t=25,b=5))
      
      p3 <- add_heatmap(p3, z=aval3[[i]]$z, x=aval3[[i]]$x,  y=aval3[[i]]$y, visible = aval3[[i]]$visible, hoverinfo="none", 
                        zmax=max.C, zmin=min.C, colorbar=list(title = "Right", x = 1, xanchor="left", y = 0.5, yanchor="top", thicknessmode="fraction", thickness=0.05, lenmode="fraction", len=0.5))
      p3 <- p3 %>% layout(xaxis = ax, yaxis = ay, title=list(text=""), margin=list(l=5,t=25,b=5))
      
      step <- list(args = list('visible', rep(FALSE, length(aval1))), label=x$times[i], method = 'restyle')
      step$args[[2]][i] = TRUE  
      steps[[i]] = step 
    }  
    p <- subplot(p1,p2,p3, margin=0.02) %>%
      layout(sliders = list(list(currentvalue = list(prefix = "Time coordinate: "), steps = steps)))
    return(p)
  }
}



### Vizualization based on slabs/slices of the data ###
#######################################################

intensitySTPP.slice <- function(X, n.slices, sigma=NULL){
  # X .......... marked point pattern, marks contain the temporal coordinates of observed events
  # n.slices ... number of slabs/slices to be considered
  # sigma ...... bandwidth parameter for the spatial kernel
  
  Xspace <- unmark(X)
  
  # Determine times where to cut slices
  t <- marks(X)
  ot <- sort(marks(X))
  t.sep <- as.numeric(quantile(ot, probs=seq(from=0, to=1, length.out=(n.slices+1))))
  
  # Create individual subpatterns as spatial projections of slices
  subpatterns <- list()
  for (k in 1:n.slices){
    if (k==1){t.keep <- (t <= t.sep[k+1])} else {t.keep <- (t <= t.sep[k+1]) & (t > t.sep[k])}
    subpatterns[[k]] <- ppp(x=Xspace$x[t.keep], y=Xspace$y[t.keep], window=X$window)
  }
  
  print("Splitting times for creating individual subpatterns (slices):")
  print(t.sep)
  
  print("Number of points in individual subpatterns (slices):")
  for (k in 1:n.slices){
    cat(npoints(subpatterns[[k]])," ")
  }
  
  # Estimate intensity functions of individual subpatterns
  intensities.P <- list()
  for (k in 1:n.slices){intensities.P[[k]] <- density.ppp(subpatterns[[k]], kernel="gaussian", sigma=sigma)}
  
  # Compute mean intensity function of the individual subpatterns  
  intensity.M <- intensities.P[[1]]
  for (k in 2:n.slices){
    aux <- intensities.P[[k]]
    intensity.M <- eval.im(intensity.M + aux)
  }
  intensity.M <- eval.im(intensity.M / n.slices)
  
  # Estimate pointwise standard deviations
  intensity.SD <- eval.im(intensity.M-intensity.M)
  for (k in 1:n.slices){
    aux <- intensities.P[[k]]
    intensity.SD <- eval.im(intensity.SD + (aux-intensity.M)^2)
  }
  intensity.SD <- eval.im(sqrt(intensity.SD / (n.slices-1)))
  
  # Compute coefficient of variation
  intensity.CV <- eval.im(intensity.SD / intensity.M)
  
  # Compute standardized values of intensity functions
  intensities.std <- list()
  for (k in 1:n.slices){
    aux <- intensities.P[[k]]
    intensities.std[[k]] <- eval.im((aux-intensity.M)/intensity.SD)
  }
  
  return(list(intensities.slices=intensities.P, intensity.mean=intensity.M,
              intensity.SD=intensity.SD, intensity.CV=intensity.CV,
              intensities.std=intensities.std))
}

plot.intensitySTPP.slice <- function(x, same.range=TRUE){
  # x           ... output of function "IntensityST.slice"
  # same.range  ... logical, should the same color scale be used for all subplots?
  
  n.patterns <- length(x$intensities.slices)
  
  if ((n.patterns==1)){
    stop("Interactive plotting not available for a single subpattern, switching to static plotting instead.")
  }
  
  # Plotting, either with the same color scale for all subpatterns
  # or with the color scale fully exploited in each subpattern
  
  ax <- list(
    visible = FALSE
  )
  ay <- list(
    visible = FALSE,
    scaleanchor = "x"
  )
  
  # Estimated intensity functions of individual slices
  aval1 <- list()
  for (k in 1:n.patterns){
    aval1[[k]] <- list(visible = FALSE, z=x$intensities.slices[[k]]$v, x=x$intensities.slices[[k]]$xcol, y=x$intensities.slices[[k]]$yrow)
  }
  aval1[[1]]$visible = TRUE
  
  # Mean of estimated intensity functions of individual slices
  aval2 <- list()
  for (k in 1:n.patterns){
    aval2[[k]] <- list(visible = FALSE, z=x$intensity.mean$v, x=x$intensity.mean$xcol, y=x$intensity.mean$yrow)
  }
  aval2[[1]]$visible = TRUE
  
  # Estimated intensity functions of individual slices
  aval3 <- list()
  for (k in 1:n.patterns){
    aval3[[k]] <- list(visible = FALSE, z=x$intensities.std[[k]]$v, x=x$intensities.std[[k]]$xcol, y=x$intensities.std[[k]]$yrow)
  }
  aval3[[1]]$visible = TRUE
  
  
  if (!same.range){
    steps <- list()
    p <- p1 <- p2 <- p3 <- plot_ly()
    for (i in 1:n.patterns) {
      zmax <- max(max(aval1[[i]]$z, na.rm=TRUE),max(aval2[[i]]$z, na.rm=TRUE), na.rm=TRUE)
      zmin <- min(min(aval1[[i]]$z, na.rm=TRUE),min(aval2[[i]]$z, na.rm=TRUE), na.rm=TRUE)
      
      p1 <- add_heatmap(p1, z=aval1[[i]]$z, x=aval1[[i]]$x,  y=aval1[[i]]$y, visible = aval1[[i]]$visible, hoverinfo="none", 
                        zmax=zmax, zmin=zmin, colorbar=list(title = "Left, middle", x = 1, xanchor="left", y = 1, yanchor="top", thicknessmode="fraction", thickness=0.05, lenmode="fraction", len=0.5))
      p1 <- p1 %>% layout(xaxis = ax, yaxis = ay, title=list(text=""), margin=list(l=5,t=25,b=5))
      
      p2 <- add_heatmap(p2, z=aval2[[i]]$z, x=aval2[[i]]$x,  y=aval2[[i]]$y, visible = aval2[[i]]$visible, hoverinfo="none",
                        zmax=zmax, zmin=zmin, colorbar=list(title = "Left, middle", x = 1, xanchor="left", y = 1, yanchor="top", thicknessmode="fraction", thickness=0.05, lenmode="fraction", len=0.5))
      p2 <- p2 %>% layout(xaxis = ax, yaxis = ay, title=list(text=""), margin=list(l=5,t=25,b=5))
      
      p3 <- add_heatmap(p3, z=aval3[[i]]$z, x=aval3[[i]]$x,  y=aval3[[i]]$y, visible = aval3[[i]]$visible, hoverinfo="none", 
                        colorbar=list(title = "Right", x = 1, xanchor="left", y = 0.5, yanchor="top", thicknessmode="fraction", thickness=0.05, lenmode="fraction", len=0.5))
      p3 <- p3 %>% layout(xaxis = ax, yaxis = ay, title=list(text=""), margin=list(l=5,t=25,b=5))
      
      step <- list(args = list('visible', rep(FALSE, length(aval1))), label=i, method = 'restyle')
      step$args[[2]][i] = TRUE  
      steps[[i]] = step 
    }  
    p <- subplot(p1,p2,p3, margin=0.02) %>%
      layout(sliders = list(list(currentvalue = list(prefix = "Slice: "), steps = steps)))
    return(p)
  }
  
  if (same.range){
    max.valsA <- max.valsB <- max.valsC <- rep(NA, times=n.patterns)
    min.valsA <- min.valsB <- min.valsC <- rep(NA, times=n.patterns)
    for (k in 1:n.patterns){
      max.valsA[k] <- max(aval1[[k]]$z, na.rm=TRUE)
      max.valsB[k] <- max(aval2[[k]]$z, na.rm=TRUE)
      max.valsC[k] <- max(aval3[[k]]$z, na.rm=TRUE)
      min.valsA[k] <- min(aval1[[k]]$z, na.rm=TRUE)
      min.valsB[k] <- min(aval2[[k]]$z, na.rm=TRUE)
      min.valsC[k] <- min(aval3[[k]]$z, na.rm=TRUE)
    }
    max.AB <- max(c(max.valsA,max.valsB))
    min.AB <- min(c(min.valsA,min.valsB))
    max.C <- max(max.valsC)
    min.C <- min(min.valsC)
    
    steps <- list()
    p <- p1 <- p2 <- p3 <- plot_ly()
    for (i in 1:n.patterns) {
      p1 <- add_heatmap(p1, z=aval1[[i]]$z, x=aval1[[i]]$x,  y=aval1[[i]]$y, visible = aval1[[i]]$visible, hoverinfo="none", 
                        zmax=max.AB, zmin=min.AB, colorbar=list(title = "Left, middle", x = 1, xanchor="left", y = 1, yanchor="top", thicknessmode="fraction", thickness=0.05, lenmode="fraction", len=0.5))
      p1 <- p1 %>% layout(xaxis = ax, yaxis = ay, title=list(text=""), margin=list(l=5,t=25,b=5))
      
      p2 <- add_heatmap(p2, z=aval2[[i]]$z, x=aval2[[i]]$x,  y=aval2[[i]]$y, visible = aval2[[i]]$visible, hoverinfo="none",
                        zmax=max.AB, zmin=min.AB, colorbar=list(title = "Left, middle", x = 1, xanchor="left", y = 1, yanchor="top", thicknessmode="fraction", thickness=0.05, lenmode="fraction", len=0.5))
      p2 <- p2 %>% layout(xaxis = ax, yaxis = ay, title=list(text=""), margin=list(l=5,t=25,b=5))
      
      p3 <- add_heatmap(p3, z=aval3[[i]]$z, x=aval3[[i]]$x,  y=aval3[[i]]$y, visible = aval3[[i]]$visible, hoverinfo="none", 
                        zmax=max.C, zmin=min.C, colorbar=list(title = "Right", x = 1, xanchor="left", y = 0.5, yanchor="top", thicknessmode="fraction", thickness=0.05, lenmode="fraction", len=0.5))
      p3 <- p3 %>% layout(xaxis = ax, yaxis = ay, title=list(text=""), margin=list(l=5,t=25,b=5))
      
      step <- list(args = list('visible', rep(FALSE, length(aval1))), label=i, method = 'restyle')
      step$args[[2]][i] = TRUE  
      steps[[i]] = step 
    }  
    p <- subplot(p1,p2,p3, margin=0.02) %>%
      layout(sliders = list(list(currentvalue = list(prefix = "Slice: "), steps = steps)))
    return(p)
  }
}
### Estimation of summary statistics for space-time point patterns (STPPs),
### to be used for stochastic reconstruction of STPPs.

### Currently the estimation is implemented for:
###   - space-time intensity function in a pixel grid (kernel estimation in separable or non-separable form)
###   - space-time intensity function in the observed points (kernel estimation in separable or non-separable form)
###   - inhomogeneous space-time K-function
###   - inhomogeneous space-time L-function
###   - raw versions of the space-time D_k-functions, ignoring inhomogeneity and edge-effects
###     (only relevant for the stochastic reconstruction!), more details are given below.

### Implemented by JiĹ™Ă­ DvoĹ™Ăˇk

### Version 28.3.2021

### Codes accompany the paper "Testing the first-order separability hypothesis for spatio-temporal point patterns"
### by M. Ghorbani, N. Vafaei, J. DvoĹ™Ăˇk and M. MyllymĂ¤ki,
### published in the Computational Statistics and Data Analysis journal 161, 107245.
### https://doi.org/10.1016/j.csda.2021.107245



library(spatstat)


########################################################################
### Temporal kernels for estimation of space-time intensity function ###
########################################################################

# Gaussian kernel
ker.time <- function(si, t, sigma.t){
  return(dnorm(x=si, mean=t, sd=sigma.t))
}

# Epanechnikov kernel
# ker.time <- function(si, t, bw){
#   aux <- pmin(abs(si - t)/bw,1)
#   return(15/(16*bw)*(1-aux^2)^2)
# }



#####################################
### Space-time intensity function ###
#####################################

### Estimate the space-time intensity function (NON-SEPARABLE) in pixel grid

intensityST.nonsep.grid <- function(X, sigma.t=NULL, sigma.s=NULL, times=NULL, start.time=NULL, end.time=NULL){
  # X .......... marked point pattern, marks contain the temporal coordinates of observed events
  # sigma.t .... bandwidth for the temporal kernel
  # sigma.s .... bandwidth parameter for the spatial kernel
  # times ...... vector of values t_i for which to compute \hat{\lambda}(.,t_i)
  # start.time, end.time ... endpoints of T, used for computing temporal edge-correction factors
  
  if (!is.numeric(times)){return(warning("times must be a numeric vector."))}
  
  if (is.null(start.time)){start.time <- min(X$marks)}
  if (is.null(end.time)){end.time <- max(X$marks)}
  
  # Spatial projection process
  Xspace <- unmark(X)
  
  # Output, empty so far
  est3d <- list()
  
  # Edge-correction factors, temporal part
  edgewt.t <- pnorm((end.time-X$marks)/sigma.t) - pnorm((start.time-X$marks)/sigma.t)
  # use 1/edgewt.t !
  
  for (i in 1:length(times)){
    t1 <- times[i]
    
    # Weights for the kernel estimation in 2D
    weights1 <- ker.time(si=X$marks, t=t1, sigma.t=sigma.t)/edgewt.t
    
    # Kernel estimation in 2D with weights
    dens <- density.ppp(Xspace, kernel="gaussian", weights=weights1, sigma=sigma.s, diggle=TRUE)
    est3d[[i]] <- dens
  }
  
  return(est3d)
}



### Estimate the space-time intensity function (SEPARABLE) in pixel grid

intensityST.sep.grid <- function(X, sigma.t=NULL, sigma.s=NULL, times=NULL, start.time=NULL, end.time=NULL){
  # X .......... marked point pattern, marks contain the temporal coordinates of observed events
  # sigma.t .... bandwidth for the temporal kernel
  # sigma.s .... bandwidth parameter for the spatial kernel
  # times ...... vector of values t_i for which to compute \hat{\lambda}(.,t_i)
  # start.time, end.time ... endpoints of T, used for computing temporal edge-correction factors
  
  if (!is.numeric(times)){return(warning("times must be a numeric vector."))}
  
  if (is.null(start.time)){start.time <- min(X$marks)}
  if (is.null(end.time)){end.time <- max(X$marks)}
  
  # Spatial projection process
  Xspace <- unmark(X)
  
  # Intensity function of the spatial projection process
  int.space <- density.ppp(Xspace, kernel="gaussian", sigma=sigma.s, diggle=TRUE)
  
  # Edge-correction factors, temporal part
  edgewt.t <- pnorm((end.time-X$marks)/sigma.t) - pnorm((start.time-X$marks)/sigma.t)
  # use 1/edgewt.t !
  
  # Intensity function of the temporal projection process
  int.time <- rep(NA, times=length(times))
  for (i in 1:length(times)){
    t1 <- times[i]
    weights1 <- ker.time(si=X$marks, t=t1, sigma.t=sigma.t)/edgewt.t
    int.time[i] <- sum(weights1)
  }
  
  # Output, empty so far
  est.sep <- list()
  
  # Take product of spatial and temporal estimates
  for (i in 1:length(times)){est.sep[[i]] <- eval.im(int.space*int.time[i]/X$n)}
  
  return(est.sep)
}



### Estimate space-time intensity function (NON-SEPARABLE) at points

intensityST.nonsep.points <- function(X, sigma.t=NULL, sigma.s=NULL, start.time=NULL, end.time=NULL){
  # X .......... marked point pattern, marks contain the temporal coordinates of observed events
  # sigma.t .... bandwidth for the temporal kernel
  # sigma.s .... bandwidth parameter for the spatial kernel
  # start.time, end.time ... endpoints of T, used for computing temporal edge-correction factors
  
  if (is.null(start.time)){start.time <- min(X$marks)}
  if (is.null(end.time)){end.time <- max(X$marks)}
  
  # Spatial projection process
  Xspace <- unmark(X)
  
  # Event times (temporal projection process)
  ts <- marks(X)
  
  
  # Edge-correction factors, spatial part
  win <- Xspace$window
  if (win$type == "rectangle") { # rectangular window
    xr <- win$xrange
    yr <- win$yrange
    xx <- Xspace$x
    yy <- Xspace$y
    xprob <- pnorm(xr[2], mean = xx, sd = sigma.s) - pnorm(xr[1], mean = xx, sd = sigma.s)
    yprob <- pnorm(yr[2], mean = yy, sd = sigma.s) - pnorm(yr[1], mean = yy, sd = sigma.s)
    edgewt.s <- xprob * yprob
  } else { # non-rectangular window
    edg <- second.moment.calc(Xspace, sigma = sigma.s, kernel = "gaussian", scalekernel = TRUE, what = "edge")
    edgewt.s <- safelookup(edg, Xspace, warn = FALSE)
  }
  # use 1/edgewt.s !
  
  
  # Edge-correction factors, temporal part
  edgewt.t <- pnorm((end.time-ts)/sigma.t) - pnorm((start.time-ts)/sigma.t)
  # use 1/edgewt.t !
  
  
  # Estimation of the intensity function
  int.out <- rep(NA, times=X$n)
  
  for (i in 1:X$n){
    xi <- X$x[i]
    yi <- X$y[i]
    ti <- X$marks[i]
    
    k1 <- dnorm(X$x, mean=xi, sd=sigma.s)*dnorm(X$y, mean=yi, sd=sigma.s)/edgewt.s
    k2 <- dnorm(X$marks, mean=ti, sd=sigma.t)/edgewt.t
    
    int.out[i] <- sum(k1*k2)
  }
  
  return(int.out)
}



### Estimate space-time intensity function (SEPARABLE) at points

intensityST.sep.points <- function(X, sigma.t=NULL, sigma.s=NULL, start.time=NULL, end.time=NULL){
  # X .......... marked point pattern, marks contain the temporal coordinates of observed events
  # sigma.t .... bandwidth for the temporal kernel
  # sigma.s .... bandwidth parameter for the spatial kernel
  # start.time, end.time ... endpoints of T, used for computing temporal edge-correction factors
  
  if (is.null(start.time)){start.time <- min(X$marks)}
  if (is.null(end.time)){end.time <- max(X$marks)}
  
  # Spatial projection process
  Xspace <- unmark(X)
  
  # Event times (temporal projection process)
  ts <- marks(X)
  
  
  # Edge-correction factors, spatial part
  win <- Xspace$window
  if (win$type == "rectangle") { # rectangular window
    xr <- win$xrange
    yr <- win$yrange
    xx <- Xspace$x
    yy <- Xspace$y
    xprob <- pnorm(xr[2], mean = xx, sd = sigma.s) - pnorm(xr[1], mean = xx, sd = sigma.s)
    yprob <- pnorm(yr[2], mean = yy, sd = sigma.s) - pnorm(yr[1], mean = yy, sd = sigma.s)
    edgewt.s <- xprob * yprob
  } else { # non-rectangular window
    edg <- second.moment.calc(Xspace, sigma = sigma.s, kernel = "gaussian", scalekernel = TRUE, what = "edge")
    edgewt.s <- safelookup(edg, Xspace, warn = FALSE)
  }
  # use 1/edgewt.s !
  
  
  # Edge-correction factors, temporal part
  edgewt.t <- pnorm((end.time-ts)/sigma.t) - pnorm((start.time-ts)/sigma.t)
  # use 1/edgewt.t !
  
  
  # Estimation of the intensity function
  int.out <- rep(NA, times=X$n)
  
  for (i in 1:X$n){
    xi <- X$x[i]
    yi <- X$y[i]
    ti <- X$marks[i]
    
    k1 <- dnorm(X$x, mean=xi, sd=sigma.s)*dnorm(X$y, mean=yi, sd=sigma.s)/edgewt.s
    k2 <- dnorm(X$marks, mean=ti, sd=sigma.t)/edgewt.t
    
    int.out[i] <- sum(k1)*sum(k2)/X$n
  }
  
  return(int.out)
}



############################################
### Space-time K-function and L-function ###
############################################

# The estimation procedure incorporates the renormalization term as in the spatstat "Kinhom" function

Krt <- function (X, lambda, Rmax=NULL, Tmax=NULL, grid=16, start.time=NULL, end.time=NULL, gW=NULL){
  # X .......... marked point pattern, marks contain the temporal coordinates of observed events
  # Rmax ....... upper limit for the K(r,t) function in the first argument
  # Tmax ....... upper limit for the K(r,t) function in the second argument
  # grid ....... number of grid points in each coordinate of K(r,t)
  # start.time, end.time ... endpoints of T
  # gW ......... set covariance of W (if already computed)
  
  if ((!is.numeric(Tmax)) || (Tmax <= 0)) stop("positive value Tmax must be provided!")
  if ((!is.numeric(Rmax)) || (Rmax <= 0)) stop("positive value Rmax must be provided!")
  if (is.null(start.time)){start.time <- min(X$marks)}
  if (is.null(end.time)){end.time <- max(X$marks)}
  
  # Vector of arguments "r" at which to compute the estimates
  r <- seq(from=Rmax/grid, to=Rmax, length.out=grid)
  
  # Renormalization factor as in the spatstat function Kinhom, with normpower=2
  renorm <- (area(X$window)*(end.time-start.time)/sum(1/lambda))^2
  # print(sum(1/lambda))
  # print(renorm)
  # renorm <- 1
  
  reciplambda <- 1/lambda
  
  # Spatial projection process
  Xspace <- unmark(X)
  
  # Determine the pairs of points which are spatially close
  closes <- closepairs(Xspace,Rmax)
  I <- closes$i
  J <- closes$j
  
  # Determine the time difference in all the pairs
  Xt <- as.numeric(X$marks)
  XtI <- Xt[I]
  XtJ <- Xt[J]
  dIJtime <- abs(XtI - XtJ)
  
  # Matrix to store the estimates
  output <- matrix(NA, ncol=grid, nrow=grid)
  
  # Compute the set covariance of W, if not already computed
  if (is.null(gW)) gW <- setcov(Window(X))
  
  # Estimate K(r,t_i) for t_i = Tmax*i/grid
  for (i in 1:grid){
    # for the moment keep only the pairs of points that are also close in the temporal domain
    keep <- (dIJtime <= Tmax*i/grid)
    
    I1   <- closes$i[keep]
    J1   <- closes$j[keep]
    dIJ1 <- closes$d[keep]
    dIJtime1 <- dIJtime[keep]
    XtI1 <- Xt[I1]
    XtJ1 <- Xt[J1]
    
    # Estimate the values using the weighted histogram approach as in the Kinhom function from spatstat
    YI <- Xspace[I1]
    YJ <- Xspace[J1]
    wIJ <- reciplambda[I1]*reciplambda[J1]
    edgewts <- edge.Trans(YI,YJ,paired=T, gW=gW)
    edgewtt <- 1 + ((XtI1 + dIJtime1 > end.time)|(XtI1 - dIJtime1 < start.time))
    allweight <- edgewts*edgewtt*wIJ
    
    wh <- whist(dIJ1,c(0,r),allweight)
    Kspac <- cumsum(wh)/(area(X$window)*(end.time-start.time))
    
    output[i,] <- Kspac
    # t varies across rows,
    # r varies across columns
  }
  
  # Apply the renormalization and return the estimates
  return( renorm*output )
}


# Space-time L-function as the square root of the space-time K-function

Lrt <- function (X, lambda, Tmax=NULL, Rmax=NULL, grid=16, start.time=NULL, end.time=NULL, gW=NULL){
  # X .......... marked point pattern, marks contain the temporal coordinates of observed events
  # Rmax ... upper limit for the L(r,t) function in the first argument
  # Tmax ... upper limit for the L(r,t) function in the second argument
  # grid ... number of grid points in each coordinate of L(r,t)
  # start.time, end.time ... endpoints of T
  # gW ......... set covariance of W (if already computed)
  
  # Estimate the space-time K-function
  aux <- Krt(X=X, lambda=lambda, Tmax=Tmax, Rmax=Rmax, grid=grid, start.time=start.time, end.time=end.time, gW=gW)
  
  # Take square root of the space-time K-function
  return(sqrt(aux))
}



################################
### Space-time D_k-functions ###
################################

# D_k(r,t) is the raw estimate of probability that a point of the process has
# at least k neighbours within spatial distance r and temporal distance t.

Drt <- function (X, Kmax=NULL, Tmax=NULL, Rmax=NULL, grid=16){
  # X .......... marked point pattern, marks contain the temporal coordinates of observed events
  # Rmax ... upper limit for the D_k(r,t) function in the first argument
  # Tmax ... upper limit for the D_k(r,t) function in the second argument
  # grid ... number of grid points in each coordinate of D_k(r,t)
  
  if ((!is.numeric(Kmax)) || (Kmax <= 0)) stop("positive value Kmax must be provided!")
  if ((!is.numeric(Tmax)) || (Tmax <= 0)) stop("positive value Tmax must be provided!")
  if ((!is.numeric(Rmax)) || (Rmax <= 0)) stop("positive value Rmax must be provided!")
  
  # Make sure Kmax is a positive integer
  Kmax <- ceiling(Kmax)
  
  # Vectors of arguments at which to compute the estimates
  r <- seq(from=Rmax/grid, to=Rmax, length.out=grid)
  t <- seq(from=Tmax/grid, to=Tmax, length.out=grid)
  
  # Spatial projection process
  Xspace <- unmark(X)
  
  # Determine the pairs of points which are spatially close
  closes <- closepairs(Xspace,Rmax)
  I <- closes$i
  J <- closes$j
  dIJspace <- closes$d
  
  # Determine the time difference in all the pairs
  Xt <- as.numeric(X$marks)
  XtI <- Xt[I]
  XtJ <- Xt[J]
  dIJtime <- abs(XtI - XtJ)
  
  # Array to store the estimates
  output <- array(NA, dim=c(grid,grid,Kmax))
  
  for (i.r in 1:grid){ 
    for (i.t in 1:grid){ 
      aux <- split((dIJspace <= r[i.r]) & (dIJtime <= t[i.t]), I, drop=FALSE)
      aux2 <- as.numeric(sapply(aux, sum))
      for (i.k in 1:Kmax){ 
        output[i.r,i.t,i.k] <- sum( aux2 >= i.k ) / X$n
      } 
    }
  }
  
  return(output)
}
### Stochastic reconstruction for space-time point patterns.

### Currently two versions are implemented, based on:
###   - space-time intensity function and space-time L-function,
###   - space-time intensity function, space-time L-function and space-time D_k-functions.

### Implemented by JiĹ™Ă­ DvoĹ™Ăˇk

### Version 28.3.2021

### Codes accompany the paper "Testing the first-order separability hypothesis for spatio-temporal point patterns"
### by M. Ghorbani, N. Vafaei, J. DvoĹ™Ăˇk and M. MyllymĂ¤ki,
### published in the Computational Statistics and Data Analysis journal 161, 107245.
### https://doi.org/10.1016/j.csda.2021.107245



library(spatstat)
library(truncnorm)


### Stochastic reconstruction procedure using space-time intensity function, L(r,t)

reconstruct.ST <- function(X, sigma.s, sigma.t, ts, start.time, end.time, Rmax, Tmax, grid, controls){
  # X ............ point pattern of spatial locations with time coordinates as marks
  # sigma.s ...... sd of the spatial Gaussian smoothing kernel
  # sigma.t ...... sd of the temporal Gaussian smoothing kernel
  # ts ........... vector of times at which to compute estimates of lambda(u,t)
  # start.time ... beginning of the time interval
  # end.time ..... end of the time interval
  # Rmax ......... upper limit for the L(r,t) function in the first argument
  # Tmax ......... upper limit for the L(r,t) function in the second argument
  # grid ......... number of grid points in each coordinate of L(r,t)
  # controls ..... list of control values, having the following elements:
  #                  max.it ... maximum number of iterations to be performed
  #                  max.reject ... maximum number of successive rejections of proposed patterns in a row
  #                  weight.L ... weight of the L-function term in the energy functional
  #                  weight.int ... weight of the intensity function term in the energy functional
  
  
  max.it <- controls$max.it
  max.reject <- controls$max.reject
  weight.L <- controls$weight.L
  weight.int <- controls$weight.int
  
  
  ### Define the energy functional ###
  ####################################
  
  energy <- function(lambda.target, L.target, lambda.act, L.act){
    EL <- sum((L.act - L.target)^2)
    Elambda <- 0
    for (i in 1:length(lambda.target)){
      aux.1 <- lambda.act[[i]]
      aux.2 <- lambda.target[[i]]
      aux.im <- eval.im((aux.1 - aux.2)^2)
      Elambda <- Elambda + integral.im(aux.im)
    }
    
    # Scaling so that the value does not depend on the number of arguments for the L(r,t) function
    EL <- EL / (nrow(L.target)*ncol(L.target))
    
    # Scaling by pixel area so that the value does not depend on the grid step
    Elambda <- Elambda*lambda.target[[1]]$xstep*lambda.target[[1]]$ystep
    
    # Scaling so that the value does not depend on the number of time slices used
    Elambda <- Elambda/length(lambda.target)
    
    # Total energy, the weights must be tuned up manually!
    res <- weight.L*EL + weight.int*Elambda 
    return(c(res,weight.L*EL,weight.int*Elambda))
  }
  
  
  ### Stochastic reconstruction ###
  #################################
  
  cat("Estimating summary characteristics of the input pattern.")
  cat("\n")
  
  n = X$n                # fixed number of points in the patterns
  count = 0              # actual number of iteration steps performed so far
  en <- matrix(NA, nrow=3, ncol=(max.it+1)) # en[1,] ..... vector of values of the energy functional for the series of patterns (for plotting)
  # en[2:3,] ... auxiliary vectors storing values of the individual terms in the energy functional
  
  same = 0               # how many iteration steps in a row without accepting a proposal
  
  # Target intensity function
  lambda.input <- intensityST.sep.grid(X, sigma.t=sigma.t, sigma.s=sigma.s, times=ts, start.time=start.time, end.time=end.time)
  
  # Target L-function
  gW <- setcov(Window(X))
  aux.lambda <- intensityST.nonsep.points(X, sigma.t=sigma.t, sigma.s=sigma.s, start.time=start.time, end.time=end.time)
  L.input <- Lrt(X, lambda=aux.lambda, Tmax=Tmax, Rmax=Rmax, start.time=start.time, end.time=end.time, grid=grid, gW=gW)
  
  
  cat("Generating the initial configuration.")
  cat("\n\n")
  
  # Initial configuration, spatial coordinates
  lambda.space.input <- density(X, sigma=sigma.s)
  Y0 <- rpoint(n=n, f=lambda.space.input)
  
  # Initial configuration, temporal coordinates
  aux.prob <- pnorm((end.time-X$marks)/sigma.t) - pnorm((start.time-X$marks)/sigma.t)
  aux.components <- sample(x=1:n, size=n, replace=TRUE, prob=aux.prob)
  aux.marks <- rtruncnorm(n=n, mean=X$marks[aux.components], sd=sigma.t, a=start.time, b=end.time)
  Y <- Y0 %mark% aux.marks
  
  # Summaries of the initial configuration
  lambda.act <- intensityST.sep.grid(Y, sigma.t=sigma.t, sigma.s=sigma.s, times=ts, start.time=start.time, end.time=end.time)
  aux.lambda <- intensityST.sep.points(Y, sigma.t=sigma.t, sigma.s=sigma.s, start.time=start.time, end.time=end.time)
  L.act <- Lrt(Y, lambda=aux.lambda, Tmax=Tmax, Rmax=Rmax, start.time=start.time, end.time=end.time, grid=grid, gW=gW)
  
  # Energy of the initial configuration
  en[,1] <- energy(lambda.target=lambda.input, L.target=L.input, lambda.act=lambda.act, L.act=L.act)[1]
  
  
  cat(paste0("Maximum number of iterations: ",max.it,". Maximum number of rejected proposals in a row: ",max.reject,"."))
  cat("\n")
  cat("Starting interations, the ''+'' symbol indicates accepted proposal, the  ''.'' symbol indicates rejected proposal.")
  cat("\n\n")
  
  # Iteration step
  while (count < max.it & same <= max.reject){
    
    # New proposal
    i <- sample(1:n, 1)
    y <- rpoint(n=1, f=lambda.space.input)
    Ynew <- Y
    Ynew$x[i] <- y$x
    Ynew$y[i] <- y$y
    aux.comp <- sample(x=1:n, size=1, replace=TRUE, prob=aux.prob)
    aux.mark <- rtruncnorm(n=1, mean=X$marks[aux.comp], sd=sigma.t, a=start.time, b=end.time)
    Ynew$marks[i] <- aux.mark
    
    # Summaries of the proposed pattern
    lambda.act <- intensityST.sep.grid(Ynew, sigma.t=sigma.t, sigma.s=sigma.s, times=ts, start.time=start.time, end.time=end.time)
    aux.lambda <- intensityST.sep.points(Ynew, sigma.t=sigma.t, sigma.s=sigma.s, start.time=start.time, end.time=end.time)
    L.act <- Lrt(Ynew, lambda=aux.lambda, Tmax=Tmax, Rmax=Rmax, start.time=start.time, end.time=end.time, grid=grid, gW=gW)
    
    # Energy of the proposed pattern
    ennew <- energy(lambda.target=lambda.input, L.target=L.input, lambda.act=lambda.act, L.act=L.act)[1]
    
    # Accept or reject the proposal based on the energy
    if (ennew[1] <= en[1,count + 1]){
      Y <- Ynew
      en[,count + 2] <- ennew
      same <- 0} 
    else {
      en[,count + 2] <- en[,count + 1]
      same <- same + 1
    }
    count <- count + 1
    if (same==0){
      cat("+") # proposal accepted
    } else {
      cat(".") # proposal rejected
    }
    if ((count %% 100) == 0){
      cat("(finished ")
      cat(count)
      cat(" iterations so far)")
      cat("\n")
      par(mfrow=c(1,2))
      plot(en[1,1:(count+1)], type="l", main="Linear scale", ylab="Energy", xlab="Iterations", ylim=c(0,en[1,1]))
      points(x=(1:(count+1)), y=en[2,1:(count+1)], type="l", col=2) # L
      points(x=(1:(count+1)), y=en[3,1:(count+1)], type="l", col=3) # intensity
      legend("topright", legend=c("Total","L","Int"), lwd=c(1,1,1,1), col=c(1,2,3,4))
      plot(en[1,1:(count+1)], type="l", main="Logarithmic scale", log="y", ylab="Energy", xlab="Iterations", ylim=c(1,en[1]))
      points(x=(1:(count+1)), y=en[2,1:(count+1)], type="l", col=2) # L
      points(x=(1:(count+1)), y=en[3,1:(count+1)], type="l", col=3) # intensity
      legend("bottomleft", legend=c("Total","L","Int"), lwd=c(1,1,1,1), col=c(1,2,3,4))
    }
  }
  cat("\n")
  
  # Make sure the output pattern has the same window as the input pattern
  # (relevant when polygonal window of the input is changed to binary mask during reconstruction)
  Y$window <- X$window
  
  # Determine the contribution of individual parts to the energy functional, useful for tuning the weights
  energy.final <- energy(lambda.target=lambda.input, L.target=L.input, lambda.act=lambda.act, L.act=L.act)
  EL.final <- energy.final[2]
  Elambda.final <- energy.final[3]
  
  # Create list of objects to be returned
  l = list(Output = Y, Count = count, Energy = en[,1:(count+1)], EL = EL.final, Elambda = Elambda.final)  
  return(l)
}



### Stochastic reconstruction procedure using space-time intensity function, L(r,t), D_k(r,t)

reconstruct.ST2 <- function(X, sigma.s, sigma.t, ts, start.time, end.time, Rmax.L, Tmax.L,
                            Kmax, Rmax.D, Tmax.D, grid.L, grid.D, controls){
  # X ............ point pattern of spatial locations with time coordinates as marks
  # sigma.s ...... sd of the spatial Gaussian smoothing kernel
  # sigma.t ...... sd of the temporal Gaussian smoothing kernel
  # ts ........... vector of times at which to compute estimates of lambda(u,t)
  # start.time ... beginning of the time interval
  # end.time ..... end of the time interval
  # Rmax.L ....... upper limit for the K(r,t) function in the first argument
  # Tmax.L ....... upper limit for the K(r,t) function in the second argument
  # Rmax.D ....... upper limit for the D_k(r,t) function in the first argument
  # Tmax.D ....... upper limit for the D_k(r,t) function in the second argument
  # Kmax ......... D_k-functions are considered for k = 1, ..., Kmax
  # grid.L ....... number of grid points in each coordinate of L(r,t)
  # grid.D ....... number of grid points in each coordinate of D_k(r,t)
  # controls ..... list of control values, having the following elements:
  #                  max.it ... maximum number of iterations to be performed
  #                  max.reject ... maximum number of successive rejections of proposed patterns in a row
  #                  weight.L ... weight of the L-function term in the energy functional
  #                  weight.int ... weight of the intensity function term in the energy functional
  #                  weight.D ... weight of the D_k-functions term in the energy functional
  
  
  max.it <- controls$max.it
  max.reject <- controls$max.reject
  weight.L <- controls$weight.L
  weight.int <- controls$weight.int
  weight.D <- controls$weight.D
  
  
  ### Define the energy functional ###
  ####################################
  
  energy = function(lambda.target, L.target, D.target, lambda.act, L.act, D.act){
    EL <- sum((L.act - L.target)^2)
    ED <- sum((D.act - D.target)^2)
    Elambda <- 0
    for (i in 1:length(lambda.target)){
      aux.1 <- lambda.act[[i]]
      aux.2 <- lambda.target[[i]]
      aux.im <- eval.im((aux.1 - aux.2)^2)
      Elambda <- Elambda + integral.im(aux.im)
    }
    
    # Scaling so that the value does not depend on the number of arguments for the L(r,t) function
    EL <- EL / (nrow(L.target)*ncol(L.target))
    
    # Scaling so that the value does not depend on the number of arguments for the D_k(r,t) function
    ED <- ED / prod(dim(D.target))
    
    # Scaling by pixel area so that the value does not depend on the grid step
    Elambda <- Elambda*lambda.target[[1]]$xstep*lambda.target[[1]]$ystep
    
    # Scaling so that the value does not depend on the number of time slices used
    Elambda <- Elambda/length(lambda.target)
    
    # Total energy, the weights must be tuned up manually!
    res <- weight.L*EL + weight.int*Elambda + weight.D*ED
    return(c(res, weight.L*EL, weight.int*Elambda, weight.D*ED))
  }
  
  
  ### Stochastic reconstruction ###
  #################################
  
  cat("Estimating summary characteristics of the input pattern.")
  cat("\n")
  
  n <- X$n                # fixed number of points in the patterns
  count <- 0              # actual number of iteration steps performed so far
  same <- 0               # how many iteration steps in a row without accepting a proposal
  en <- matrix(NA, nrow=4, ncol=(max.it+1)) # en[1,] ..... vector of values of the energy functional for the series of patterns (for plotting)
  # en[2:4,] ... auxiliary vectors storing values of the individual terms in the energy functional
  
  # Target intensity function
  lambda.input <- intensityST.sep.grid(X, sigma.t=sigma.t, sigma.s=sigma.s, times=ts, start.time=start.time, end.time=end.time)
  
  # Target L-function
  gW <- setcov(Window(X))
  aux.lambda <- intensityST.nonsep.points(X, sigma.t=sigma.t, sigma.s=sigma.s, start.time=start.time, end.time=end.time)
  L.input <- Lrt(X, lambda=aux.lambda, Tmax=Tmax.L, Rmax=Rmax.L, start.time=start.time, end.time=end.time, grid=grid.L, gW=gW)
  
  # Target D_k-functions
  D.input <- Drt(X, Kmax=Kmax, Tmax=Tmax.D, Rmax=Rmax.D, grid=grid.D)
  
  
  cat("Generating the initial configuration.")
  cat("\n\n")
  
  # Initial configuration, spatial coordinates
  lambda.space.input <- density(X, sigma=sigma.s)
  Y0 <- rpoint(n=n, f=lambda.space.input)
  
  # Initial configuration, temporal coordinates
  aux.prob <- pnorm((end.time-X$marks)/sigma.t) - pnorm((start.time-X$marks)/sigma.t)
  aux.components <- sample(x=1:n, size=n, replace=TRUE, prob=aux.prob)
  aux.marks <- rtruncnorm(n=n, mean=X$marks[aux.components], sd=sigma.t, a=start.time, b=end.time)
  Y <- Y0 %mark% aux.marks
  
  # Summaries of the initial configuration
  lambda.act <- intensityST.sep.grid(Y, sigma.t=sigma.t, sigma.s=sigma.s, times=ts, start.time=start.time, end.time=end.time)
  aux.lambda <- intensityST.sep.points(Y, sigma.t=sigma.t, sigma.s=sigma.s, start.time=start.time, end.time=end.time)
  L.act <- Lrt(Y, lambda=aux.lambda, Tmax=Tmax.L, Rmax=Rmax.L, start.time=start.time, end.time=end.time, grid=grid.L, gW=gW)
  D.act <- Drt(Y, Kmax=Kmax, Tmax=Tmax.D, Rmax=Rmax.D, grid=grid.D)
  
  # Energy of the initial configuration
  en[,1] <- energy(lambda.target=lambda.input, L.target=L.input, D.target=D.input, lambda.act=lambda.act, L.act=L.act, D.act=D.act)
  
  
  cat(paste0("Maximum number of iterations: ",max.it,". Maximum number of rejected proposals in a row: ",max.reject,"."))
  cat("\n")
  cat("Starting interations, the ''+'' symbol indicates accepted proposal, the  ''.'' symbol indicates rejected proposal.")
  cat("\n\n")
  
  # Iteration step
  while (count < max.it & same <= max.reject){
    
    # New proposal
    i <- sample(1:n, 1)
    y <- rpoint(n=1, f=lambda.space.input)
    Ynew <- Y
    Ynew$x[i] <- y$x
    Ynew$y[i] <- y$y
    aux.comp <- sample(x=1:n, size=1, replace=TRUE, prob=aux.prob)
    aux.mark <- rtruncnorm(n=1, mean=X$marks[aux.comp], sd=sigma.t, a=start.time, b=end.time)
    Ynew$marks[i] <- aux.mark
    
    # Summaries of the proposed pattern
    lambda.act <- intensityST.sep.grid(Ynew, sigma.t=sigma.t, sigma.s=sigma.s, times=ts, start.time=start.time, end.time=end.time)
    aux.lambda <- intensityST.sep.points(Ynew, sigma.t=sigma.t, sigma.s=sigma.s, start.time=start.time, end.time=end.time)
    L.act <- Lrt(Ynew, lambda=aux.lambda, Tmax=Tmax.L, Rmax=Rmax.L, start.time=start.time, end.time=end.time, grid=grid.L, gW=gW)
    D.act <- Drt(Ynew, Kmax=Kmax, Tmax=Tmax.D, Rmax=Rmax.D, grid=grid.D)
    
    # Energy of the proposed pattern
    ennew <- energy(lambda.target=lambda.input, L.target=L.input, D.target=D.input, lambda.act=lambda.act, L.act=L.act, D.act=D.act)
    
    # Accept or reject the proposal based on the energy
    if (ennew[1] <= en[1,count + 1]){
      Y <- Ynew
      en[,count + 2] <- ennew
      same <- 0} 
    else {
      en[,count + 2] <- en[,count + 1]
      same <- same + 1
    }
    count <- count + 1
    if (same==0){
      cat("+") # proposal accepted
    } else {
      cat(".") # proposal rejected
    }
    if ((count %% 100) == 0){
      cat("(finished ")
      cat(count)
      cat(" iterations so far)")
      cat("\n")
      par(mfrow=c(1,2))
      plot(en[1,1:(count+1)], type="l", main="Linear scale", ylab="Energy", xlab="Iterations", ylim=c(0,en[1,1]))
      points(x=(1:(count+1)), y=en[2,1:(count+1)], type="l", col=2) # L
      points(x=(1:(count+1)), y=en[3,1:(count+1)], type="l", col=3) # intensity
      points(x=(1:(count+1)), y=en[4,1:(count+1)], type="l", col=4) # D_k
      legend("topright", legend=c("Total","L","Int","D_k"), lwd=c(1,1,1,1), col=c(1,2,3,4))
      plot(en[1,1:(count+1)], type="l", main="Logarithmic scale", log="y", ylab="Energy", xlab="Iterations", ylim=c(min(en[,count+1]),en[1,1]))
      points(x=(1:(count+1)), y=en[2,1:(count+1)], type="l", col=2) # L
      points(x=(1:(count+1)), y=en[3,1:(count+1)], type="l", col=3) # intensity
      points(x=(1:(count+1)), y=en[4,1:(count+1)], type="l", col=4) # D_k
      legend("bottomleft", legend=c("Total","L","Int","D_k"), lwd=c(1,1,1,1), col=c(1,2,3,4))
    }
  }
  cat("\n")
  
  # Make sure the output pattern has the same window as the input pattern
  # (relevant when polygonal window of the input is changed to binary mask during reconstruction)
  Y$window <- X$window
  
  # Determine the contribution of individual parts to the energy functional, useful for tuning the weights
  energy.final <- energy(lambda.target=lambda.input, L.target=L.input, D.target=D.input, lambda.act=lambda.act, L.act=L.act, D.act=D.act)
  EL.final <- energy.final[2]
  Elambda.final <- energy.final[3]
  ED.final <- energy.final[4]
  
  # Create list of objects to be returned
  l = list(Output = Y, Count = count, Energy = en[,1:(count+1)], EL = EL.final, Elambda = Elambda.final, ED = ED.final)  
  return(l)
}



reconstruct.ST2.batch <- function(X, sigma.s, sigma.t, ts, start.time, end.time, Rmax.L, Tmax.L,
                                  Kmax, Rmax.D, Tmax.D, grid.L, grid.D, controls){
  # X ............ point pattern of spatial locations with time coordinates as marks
  # sigma.s ...... sd of the spatial Gaussian smoothing kernel
  # sigma.t ...... sd of the temporal Gaussian smoothing kernel
  # ts ........... vector of times at which to compute estimates of lambda(u,t)
  # start.time ... beginning of the time interval
  # end.time ..... end of the time interval
  # Rmax.L ....... upper limit for the K(r,t) function in the first argument
  # Tmax.L ....... upper limit for the K(r,t) function in the second argument
  # Rmax.D ....... upper limit for the D_k(r,t) function in the first argument
  # Tmax.D ....... upper limit for the D_k(r,t) function in the second argument
  # Kmax ......... D_k-functions are considered for k = 1, ..., Kmax
  # grid.L ....... number of grid points in each coordinate of L(r,t)
  # grid.D ....... number of grid points in each coordinate of D_k(r,t)
  # controls ..... list of control values, having the following elements:
  #                  max.it ... maximum number of iterations to be performed
  #                  max.reject ... maximum number of successive rejections of proposed patterns in a row
  #                  weight.L ... weight of the L-function term in the energy functional
  #                  weight.int ... weight of the intensity function term in the energy functional
  #                  weight.D ... weight of the D_k-functions term in the energy functional
  
  
  max.it <- controls$max.it
  max.reject <- controls$max.reject
  weight.L <- controls$weight.L
  weight.int <- controls$weight.int
  weight.D <- controls$weight.D
  
  
  ### Define the energy functional ###
  ####################################
  
  energy = function(lambda.target, L.target, D.target, lambda.act, L.act, D.act){
    EL <- sum((L.act - L.target)^2)
    ED <- sum((D.act - D.target)^2)
    Elambda <- 0
    for (i in 1:length(lambda.target)){
      aux.1 <- lambda.act[[i]]
      aux.2 <- lambda.target[[i]]
      aux.im <- eval.im((aux.1 - aux.2)^2)
      Elambda <- Elambda + integral.im(aux.im)
    }
    
    # Scaling so that the value does not depend on the number of arguments for the L(r,t) function
    EL <- EL / (nrow(L.target)*ncol(L.target))
    
    # Scaling so that the value does not depend on the number of arguments for the D_k(r,t) function
    ED <- ED / prod(dim(D.target))
    
    # Scaling by pixel area so that the value does not depend on the grid step
    Elambda <- Elambda*lambda.target[[1]]$xstep*lambda.target[[1]]$ystep
    
    # Scaling so that the value does not depend on the number of time slices used
    Elambda <- Elambda/length(lambda.target)
    
    # Total energy, the weights must be tuned up manually!
    res <- weight.L*EL + weight.int*Elambda + weight.D*ED
    return(c(res, weight.L*EL, weight.int*Elambda, weight.D*ED))
  }
  
  
  ### Stochastic reconstruction ###
  #################################
  
  n <- X$n                # fixed number of points in the patterns
  count <- 0              # actual number of iteration steps performed so far
  same <- 0               # how many iteration steps in a row without accepting a proposal
  # en = numeric(max.it+1) # vector of values of the energy functional for the series of patterns (for plotting)
  en <- matrix(NA, nrow=4, ncol=(max.it+1)) # en[1,] ... vector of values of the energy functional for the series of patterns (for plotting)
  # en[2:4,] ... auxiliary vectors storing values of the individual terms in the energy functional
  
  # Target intensity function
  lambda.input <- intensityST.sep.grid(X, sigma.t=sigma.t, sigma.s=sigma.s, times=ts, start.time=start.time, end.time=end.time)
  
  # Target L-function
  gW <- setcov(Window(X))
  aux.lambda <- intensityST.nonsep.points(X, sigma.t=sigma.t, sigma.s=sigma.s, start.time=start.time, end.time=end.time)
  L.input <- Lrt(X, lambda=aux.lambda, Tmax=Tmax.L, Rmax=Rmax.L, start.time=start.time, end.time=end.time, grid=grid.L, gW=gW)
  
  # Target D_k-functions
  D.input <- Drt(X, Kmax=Kmax, Tmax=Tmax.D, Rmax=Rmax.D, grid=grid.D)
  
  
  # Initial configuration, spatial coordinates
  lambda.space.input <- density(X, sigma=sigma.s)
  Y0 <- rpoint(n=n, f=lambda.space.input)
  
  # Initial configuration, temporal coordinates
  aux.prob <- pnorm((end.time-X$marks)/sigma.t) - pnorm((start.time-X$marks)/sigma.t)
  aux.components <- sample(x=1:n, size=n, replace=TRUE, prob=aux.prob)
  aux.marks <- rtruncnorm(n=n, mean=X$marks[aux.components], sd=sigma.t, a=start.time, b=end.time)
  Y <- Y0 %mark% aux.marks
  
  # Summaries of the initial configuration
  lambda.act <- intensityST.sep.grid(Y, sigma.t=sigma.t, sigma.s=sigma.s, times=ts, start.time=start.time, end.time=end.time)
  aux.lambda <- intensityST.sep.points(Y, sigma.t=sigma.t, sigma.s=sigma.s, start.time=start.time, end.time=end.time)
  L.act <- Lrt(Y, lambda=aux.lambda, Tmax=Tmax.L, Rmax=Rmax.L, start.time=start.time, end.time=end.time, grid=grid.L, gW=gW)
  D.act <- Drt(Y, Kmax=Kmax, Tmax=Tmax.D, Rmax=Rmax.D, grid=grid.D)
  
  # Energy of the initial configuration
  en[,1] <- energy(lambda.target=lambda.input, L.target=L.input, D.target=D.input, lambda.act=lambda.act, L.act=L.act, D.act=D.act)
  
  
  # Iteration step
  while (count < max.it & same <= max.reject){
    
    # New proposal
    i <- sample(1:n, 1)
    y <- rpoint(n=1, f=lambda.space.input)
    Ynew <- Y
    Ynew$x[i] <- y$x
    Ynew$y[i] <- y$y
    aux.comp <- sample(x=1:n, size=1, replace=TRUE, prob=aux.prob)
    aux.mark <- rtruncnorm(n=1, mean=X$marks[aux.comp], sd=sigma.t, a=start.time, b=end.time)
    Ynew$marks[i] <- aux.mark
    
    # Summaries of the proposed pattern
    lambda.act <- intensityST.sep.grid(Ynew, sigma.t=sigma.t, sigma.s=sigma.s, times=ts, start.time=start.time, end.time=end.time)
    aux.lambda <- intensityST.sep.points(Ynew, sigma.t=sigma.t, sigma.s=sigma.s, start.time=start.time, end.time=end.time)
    L.act <- Lrt(Ynew, lambda=aux.lambda, Tmax=Tmax.L, Rmax=Rmax.L, start.time=start.time, end.time=end.time, grid=grid.L, gW=gW)
    D.act <- Drt(Ynew, Kmax=Kmax, Tmax=Tmax.D, Rmax=Rmax.D, grid=grid.D)
    
    # Energy of the proposed pattern
    ennew <- energy(lambda.target=lambda.input, L.target=L.input, D.target=D.input, lambda.act=lambda.act, L.act=L.act, D.act=D.act)
    
    # Accept or reject the proposal based on the energy
    if (ennew[1] <= en[1,count + 1]){
      Y <- Ynew
      en[,count + 2] <- ennew
      same <- 0} 
    else {
      en[,count + 2] <- en[,count + 1]
      same <- same + 1
    }
    count <- count + 1
  }
  
  # Make sure the output pattern has the same window as the input pattern
  # (relevant when polygonal window of the input is changed to binary mask during reconstruction)
  Y$window <- X$window
  
  # Determine the contribution of individual parts to the energy functional, useful for tuning the weights
  energy.final <- energy(lambda.target=lambda.input, L.target=L.input, D.target=D.input, lambda.act=lambda.act, L.act=L.act, D.act=D.act)
  EL.final <- energy.final[2]
  Elambda.final <- energy.final[3]
  ED.final <- energy.final[4]
  
  # Create list of objects to be returned
  l = list(Output = Y, Count = count, Energy = en[,1:(count+1)], EL = EL.final, Elambda = Elambda.final, ED = ED.final)  
  return(l)
}
### Estimation of summary statistics for space-time point patterns (STPPs),
### to be used for stochastic reconstruction of STPPs.

### Currently the estimation is implemented for:
###   - space-time intensity function in a pixel grid (kernel estimation in separable or non-separable form)
###   - space-time intensity function in the observed points (kernel estimation in separable or non-separable form)
###   - inhomogeneous space-time K-function
###   - inhomogeneous space-time L-function
###   - raw versions of the space-time D_k-functions, ignoring inhomogeneity and edge-effects
###     (only relevant for the stochastic reconstruction!), more details are given below.

### Implemented by JiĹ™Ă DvoĹ™Ăˇk

### Version 28.3.2021

### Codes accompany the paper "Testing the first-order separability hypothesis for spatio-temporal point patterns"
### by M. Ghorbani, N. Vafaei, J. DvoĹ™Ăˇk and M. MyllymĂ¤ki,
### published in the Computational Statistics and Data Analysis journal 161, 107245.
### https://doi.org/10.1016/j.csda.2021.107245



library(spatstat)


########################################################################
### Temporal kernels for estimation of space-time intensity function ###
########################################################################

# Gaussian kernel
ker.time <- function(si, t, sigma.t){
  return(dnorm(x=si, mean=t, sd=sigma.t))
}

# Epanechnikov kernel
# ker.time <- function(si, t, bw){
#   aux <- pmin(abs(si - t)/bw,1)
#   return(15/(16*bw)*(1-aux^2)^2)
# }



#####################################
### Space-time intensity function ###
#####################################

### Estimate the space-time intensity function (NON-SEPARABLE) in pixel grid

intensityST.nonsep.grid <- function(X, sigma.t=NULL, sigma.s=NULL, times=NULL, start.time=NULL, end.time=NULL){
  # X .......... marked point pattern, marks contain the temporal coordinates of observed events
  # sigma.t .... bandwidth for the temporal kernel
  # sigma.s .... bandwidth parameter for the spatial kernel
  # times ...... vector of values t_i for which to compute \hat{\lambda}(.,t_i)
  # start.time, end.time ... endpoints of T, used for computing temporal edge-correction factors
  
  if (!is.numeric(times)){return(warning("times must be a numeric vector."))}
  
  if (is.null(start.time)){start.time <- min(X$marks)}
  if (is.null(end.time)){end.time <- max(X$marks)}
  
  # Spatial projection process
  Xspace <- unmark(X)
  
  # Output, empty so far
  est3d <- list()
  
  # Edge-correction factors, temporal part
  edgewt.t <- pnorm((end.time-X$marks)/sigma.t) - pnorm((start.time-X$marks)/sigma.t)
  # use 1/edgewt.t !
  
  for (i in 1:length(times)){
    t1 <- times[i]
    
    # Weights for the kernel estimation in 2D
    weights1 <- ker.time(si=X$marks, t=t1, sigma.t=sigma.t)/edgewt.t
    
    # Kernel estimation in 2D with weights
    dens <- density.ppp(Xspace, kernel="gaussian", weights=weights1, sigma=sigma.s, diggle=TRUE)
    est3d[[i]] <- dens
  }
  
  return(est3d)
}



### Estimate the space-time intensity function (SEPARABLE) in pixel grid

intensityST.sep.grid <- function(X, sigma.t=NULL, sigma.s=NULL, times=NULL, start.time=NULL, end.time=NULL){
  # X .......... marked point pattern, marks contain the temporal coordinates of observed events
  # sigma.t .... bandwidth for the temporal kernel
  # sigma.s .... bandwidth parameter for the spatial kernel
  # times ...... vector of values t_i for which to compute \hat{\lambda}(.,t_i)
  # start.time, end.time ... endpoints of T, used for computing temporal edge-correction factors
  
  if (!is.numeric(times)){return(warning("times must be a numeric vector."))}
  
  if (is.null(start.time)){start.time <- min(X$marks)}
  if (is.null(end.time)){end.time <- max(X$marks)}
  
  # Spatial projection process
  Xspace <- unmark(X)
  
  # Intensity function of the spatial projection process
  int.space <- density.ppp(Xspace, kernel="gaussian", sigma=sigma.s, diggle=TRUE)
  
  # Edge-correction factors, temporal part
  edgewt.t <- pnorm((end.time-X$marks)/sigma.t) - pnorm((start.time-X$marks)/sigma.t)
  # use 1/edgewt.t !
  
  # Intensity function of the temporal projection process
  int.time <- rep(NA, times=length(times))
  for (i in 1:length(times)){
    t1 <- times[i]
    weights1 <- ker.time(si=X$marks, t=t1, sigma.t=sigma.t)/edgewt.t
    int.time[i] <- sum(weights1)
  }
  
  # Output, empty so far
  est.sep <- list()
  
  # Take product of spatial and temporal estimates
  for (i in 1:length(times)){est.sep[[i]] <- eval.im(int.space*int.time[i]/X$n)}
  
  return(est.sep)
}



### Estimate space-time intensity function (NON-SEPARABLE) at points

intensityST.nonsep.points <- function(X, sigma.t=NULL, sigma.s=NULL, start.time=NULL, end.time=NULL){
  # X .......... marked point pattern, marks contain the temporal coordinates of observed events
  # sigma.t .... bandwidth for the temporal kernel
  # sigma.s .... bandwidth parameter for the spatial kernel
  # start.time, end.time ... endpoints of T, used for computing temporal edge-correction factors
  
  if (is.null(start.time)){start.time <- min(X$marks)}
  if (is.null(end.time)){end.time <- max(X$marks)}
  
  # Spatial projection process
  Xspace <- unmark(X)
  
  # Event times (temporal projection process)
  ts <- marks(X)
  
  
  # Edge-correction factors, spatial part
  win <- Xspace$window
  if (win$type == "rectangle") { # rectangular window
    xr <- win$xrange
    yr <- win$yrange
    xx <- Xspace$x
    yy <- Xspace$y
    xprob <- pnorm(xr[2], mean = xx, sd = sigma.s) - pnorm(xr[1], mean = xx, sd = sigma.s)
    yprob <- pnorm(yr[2], mean = yy, sd = sigma.s) - pnorm(yr[1], mean = yy, sd = sigma.s)
    edgewt.s <- xprob * yprob
  } else { # non-rectangular window
    edg <- second.moment.calc(Xspace, sigma = sigma.s, kernel = "gaussian", scalekernel = TRUE, what = "edge")
    edgewt.s <- safelookup(edg, Xspace, warn = FALSE)
  }
  # use 1/edgewt.s !
  
  
  # Edge-correction factors, temporal part
  edgewt.t <- pnorm((end.time-ts)/sigma.t) - pnorm((start.time-ts)/sigma.t)
  # use 1/edgewt.t !
  
  
  # Estimation of the intensity function
  int.out <- rep(NA, times=X$n)
  
  for (i in 1:X$n){
    xi <- X$x[i]
    yi <- X$y[i]
    ti <- X$marks[i]
    
    k1 <- dnorm(X$x, mean=xi, sd=sigma.s)*dnorm(X$y, mean=yi, sd=sigma.s)/edgewt.s
    k2 <- dnorm(X$marks, mean=ti, sd=sigma.t)/edgewt.t
    
    int.out[i] <- sum(k1*k2)
  }
  
  return(int.out)
}



### Estimate space-time intensity function (SEPARABLE) at points

intensityST.sep.points <- function(X, sigma.t=NULL, sigma.s=NULL, start.time=NULL, end.time=NULL){
  # X .......... marked point pattern, marks contain the temporal coordinates of observed events
  # sigma.t .... bandwidth for the temporal kernel
  # sigma.s .... bandwidth parameter for the spatial kernel
  # start.time, end.time ... endpoints of T, used for computing temporal edge-correction factors
  
  if (is.null(start.time)){start.time <- min(X$marks)}
  if (is.null(end.time)){end.time <- max(X$marks)}
  
  # Spatial projection process
  Xspace <- unmark(X)
  
  # Event times (temporal projection process)
  ts <- marks(X)
  
  
  # Edge-correction factors, spatial part
  win <- Xspace$window
  if (win$type == "rectangle") { # rectangular window
    xr <- win$xrange
    yr <- win$yrange
    xx <- Xspace$x
    yy <- Xspace$y
    xprob <- pnorm(xr[2], mean = xx, sd = sigma.s) - pnorm(xr[1], mean = xx, sd = sigma.s)
    yprob <- pnorm(yr[2], mean = yy, sd = sigma.s) - pnorm(yr[1], mean = yy, sd = sigma.s)
    edgewt.s <- xprob * yprob
  } else { # non-rectangular window
    edg <- second.moment.calc(Xspace, sigma = sigma.s, kernel = "gaussian", scalekernel = TRUE, what = "edge")
    edgewt.s <- safelookup(edg, Xspace, warn = FALSE)
  }
  # use 1/edgewt.s !
  
  
  # Edge-correction factors, temporal part
  edgewt.t <- pnorm((end.time-ts)/sigma.t) - pnorm((start.time-ts)/sigma.t)
  # use 1/edgewt.t !
  
  
  # Estimation of the intensity function
  int.out <- rep(NA, times=X$n)
  
  for (i in 1:X$n){
    xi <- X$x[i]
    yi <- X$y[i]
    ti <- X$marks[i]
    
    k1 <- dnorm(X$x, mean=xi, sd=sigma.s)*dnorm(X$y, mean=yi, sd=sigma.s)/edgewt.s
    k2 <- dnorm(X$marks, mean=ti, sd=sigma.t)/edgewt.t
    
    int.out[i] <- sum(k1)*sum(k2)/X$n
  }
  
  return(int.out)
}



############################################
### Space-time K-function and L-function ###
############################################

# The estimation procedure incorporates the renormalization term as in the spatstat "Kinhom" function

Krt <- function (X, lambda, Rmax=NULL, Tmax=NULL, grid=16, start.time=NULL, end.time=NULL, gW=NULL){
  # X .......... marked point pattern, marks contain the temporal coordinates of observed events
  # Rmax ....... upper limit for the K(r,t) function in the first argument
  # Tmax ....... upper limit for the K(r,t) function in the second argument
  # grid ....... number of grid points in each coordinate of K(r,t)
  # start.time, end.time ... endpoints of T
  # gW ......... set covariance of W (if already computed)
  
  if ((!is.numeric(Tmax)) || (Tmax <= 0)) stop("positive value Tmax must be provided!")
  if ((!is.numeric(Rmax)) || (Rmax <= 0)) stop("positive value Rmax must be provided!")
  if (is.null(start.time)){start.time <- min(X$marks)}
  if (is.null(end.time)){end.time <- max(X$marks)}
  
  # Vector of arguments "r" at which to compute the estimates
  r <- seq(from=Rmax/grid, to=Rmax, length.out=grid)
  
  # Renormalization factor as in the spatstat function Kinhom, with normpower=2
  renorm <- (area(X$window)*(end.time-start.time)/sum(1/lambda))^2
  # print(sum(1/lambda))
  # print(renorm)
  # renorm <- 1
  
  reciplambda <- 1/lambda
  
  # Spatial projection process
  Xspace <- unmark(X)
  
  # Determine the pairs of points which are spatially close
  closes <- closepairs(Xspace,Rmax)
  I <- closes$i
  J <- closes$j
  
  # Determine the time difference in all the pairs
  Xt <- as.numeric(X$marks)
  XtI <- Xt[I]
  XtJ <- Xt[J]
  dIJtime <- abs(XtI - XtJ)
  
  # Matrix to store the estimates
  output <- matrix(NA, ncol=grid, nrow=grid)
  
  # Compute the set covariance of W, if not already computed
  if (is.null(gW)) gW <- setcov(Window(X))
  
  # Estimate K(r,t_i) for t_i = Tmax*i/grid
  for (i in 1:grid){
    # for the moment keep only the pairs of points that are also close in the temporal domain
    keep <- (dIJtime <= Tmax*i/grid)
    
    I1   <- closes$i[keep]
    J1   <- closes$j[keep]
    dIJ1 <- closes$d[keep]
    dIJtime1 <- dIJtime[keep]
    XtI1 <- Xt[I1]
    XtJ1 <- Xt[J1]
    
    # Estimate the values using the weighted histogram approach as in the Kinhom function from spatstat
    YI <- Xspace[I1]
    YJ <- Xspace[J1]
    wIJ <- reciplambda[I1]*reciplambda[J1]
    edgewts <- edge.Trans(YI,YJ,paired=T, gW=gW)
    edgewtt <- 1 + ((XtI1 + dIJtime1 > end.time)|(XtI1 - dIJtime1 < start.time))
    allweight <- edgewts*edgewtt*wIJ
    
    wh <- whist(dIJ1,c(0,r),allweight)
    Kspac <- cumsum(wh)/(area(X$window)*(end.time-start.time))
    
    output[i,] <- Kspac
    # t varies across rows,
    # r varies across columns
  }
  
  # Apply the renormalization and return the estimates
  return( renorm*output )
}


# Space-time L-function as the square root of the space-time K-function

Lrt <- function (X, lambda, Tmax=NULL, Rmax=NULL, grid=16, start.time=NULL, end.time=NULL, gW=NULL){
  # X .......... marked point pattern, marks contain the temporal coordinates of observed events
  # Rmax ... upper limit for the L(r,t) function in the first argument
  # Tmax ... upper limit for the L(r,t) function in the second argument
  # grid ... number of grid points in each coordinate of L(r,t)
  # start.time, end.time ... endpoints of T
  # gW ......... set covariance of W (if already computed)
  
  # Estimate the space-time K-function
  aux <- Krt(X=X, lambda=lambda, Tmax=Tmax, Rmax=Rmax, grid=grid, start.time=start.time, end.time=end.time, gW=gW)
  
  # Take square root of the space-time K-function
  return(sqrt(aux))
}



################################
### Space-time D_k-functions ###
################################

# D_k(r,t) is the raw estimate of probability that a point of the process has
# at least k neighbours within spatial distance r and temporal distance t.

Drt <- function (X, Kmax=NULL, Tmax=NULL, Rmax=NULL, grid=16){
  # X .......... marked point pattern, marks contain the temporal coordinates of observed events
  # Rmax ... upper limit for the D_k(r,t) function in the first argument
  # Tmax ... upper limit for the D_k(r,t) function in the second argument
  # grid ... number of grid points in each coordinate of D_k(r,t)
  
  if ((!is.numeric(Kmax)) || (Kmax <= 0)) stop("positive value Kmax must be provided!")
  if ((!is.numeric(Tmax)) || (Tmax <= 0)) stop("positive value Tmax must be provided!")
  if ((!is.numeric(Rmax)) || (Rmax <= 0)) stop("positive value Rmax must be provided!")
  
  # Make sure Kmax is a positive integer
  Kmax <- ceiling(Kmax)
  
  # Vectors of arguments at which to compute the estimates
  r <- seq(from=Rmax/grid, to=Rmax, length.out=grid)
  t <- seq(from=Tmax/grid, to=Tmax, length.out=grid)
  
  # Spatial projection process
  Xspace <- unmark(X)
  
  # Determine the pairs of points which are spatially close
  closes <- closepairs(Xspace,Rmax)
  I <- closes$i
  J <- closes$j
  dIJspace <- closes$d
  
  # Determine the time difference in all the pairs
  Xt <- as.numeric(X$marks)
  XtI <- Xt[I]
  XtJ <- Xt[J]
  dIJtime <- abs(XtI - XtJ)
  
  # Array to store the estimates
  output <- array(NA, dim=c(grid,grid,Kmax))
  
  for (i.r in 1:grid){ 
    for (i.t in 1:grid){ 
      aux <- split((dIJspace <= r[i.r]) & (dIJtime <= t[i.t]), I, drop=FALSE)
      aux2 <- as.numeric(sapply(aux, sum))
      for (i.k in 1:Kmax){ 
        output[i.r,i.t,i.k] <- sum( aux2 >= i.k ) / X$n
      } 
    }
  }
  
  return(output)
}



#### Funcion reconstruida 

reconstruct.ST1 <- function(X, sigma.s, sigma.t, ts, start.time, end.time, Rmax, Tmax, grid, controls){
  # X ............ point pattern of spatial locations with time coordinates as marks
  # sigma.s ...... sd of the spatial Gaussian smoothing kernel
  # sigma.t ...... sd of the temporal Gaussian smoothing kernel
  # ts ........... vector of times at which to compute estimates of lambda(u,t)
  # start.time ... beginning of the time interval
  # end.time ..... end of the time interval
  # Rmax ......... upper limit for the L(r,t) function in the first argument
  # Tmax ......... upper limit for the L(r,t) function in the second argument
  # grid ......... number of grid points in each coordinate of L(r,t)
  # controls ..... list of control values, having the following elements:
  #                  max.it ... maximum number of iterations to be performed
  #                  max.reject ... maximum number of successive rejections of proposed patterns in a row
  #                  weight.L ... weight of the L-function term in the energy functional
  #                  weight.int ... weight of the intensity function term in the energy functional
  
  
  max.it <- controls$max.it
  max.reject <- controls$max.reject
  weight.L <- controls$weight.L
  weight.int <- controls$weight.int
  
  ### Define the energy functional ###
  ####################################
  
  energy <- function(lambda.target, L.target, lambda.act, L.act){
    # Check if there are NA values in the inputs
    if (any(is.na(L.target)) || any(is.na(L.act))) {
      stop("NA values detected in L.target or L.act")
    }
    
    EL <- sum((L.act - L.target)^2)
    Elambda <- 0
    for (i in 1:length(lambda.target)){
      aux.1 <- lambda.act[[i]]
      aux.2 <- lambda.target[[i]]
      aux.im <- eval.im((aux.1 - aux.2)^2)
      Elambda <- Elambda + integral.im(aux.im)
    }
    
    # Scaling so that the value does not depend on the number of arguments for the L(r,t) function
    EL <- EL / (nrow(L.target)*ncol(L.target))
    
    # Scaling by pixel area so that the value does not depend on the grid step
    Elambda <- Elambda*lambda.target[[1]]$xstep*lambda.target[[1]]$ystep
    
    # Scaling so that the value does not depend on the number of time slices used
    Elambda <- Elambda/length(lambda.target)
    
    # Total energy, the weights must be tuned up manually!
    res <- weight.L*EL + weight.int*Elambda 
    return(c(res, weight.L*EL, weight.int*Elambda))
  }
  
  ### Stochastic reconstruction ###
  #################################
  
  cat("Estimating summary characteristics of the input pattern.\n")
  
  n = X$n                # fixed number of points in the patterns
  count = 0              # actual number of iteration steps performed so far
  en <- matrix(NA, nrow=3, ncol=(max.it+1)) # en[1,] ..... vector of values of the energy functional for the series of patterns (for plotting)
  same = 0               # how many iteration steps in a row without accepting a proposal
  
  # Target intensity function
  lambda.input <- intensityST.sep.grid(X, sigma.t=sigma.t, sigma.s=sigma.s, times=ts, start.time=start.time, end.time=end.time)
  
  # Check if lambda.input has NA values
  if (any(is.na(lambda.input))) {
    stop("NA values detected in lambda.input")
  }
  
  # Target L-function
  gW <- setcov(Window(X))
  aux.lambda <- intensityST.nonsep.points(X, sigma.t=sigma.t, sigma.s=sigma.s, start.time=start.time, end.time=end.time)
  
  # Check if aux.lambda has NA values
  if (any(is.na(aux.lambda))) {
    stop("NA values detected in aux.lambda")
  }
  
  L.input <- Lrt(X, lambda=aux.lambda, Tmax=Tmax, Rmax=Rmax, start.time=start.time, end.time=end.time, grid=grid, gW=gW)
  
  # Check if L.input has NA values
  if (any(is.na(L.input))) {
    stop("NA values detected in L.input")
  }
  
  cat("Generating the initial configuration.\n\n")
  
  # Initial configuration, spatial coordinates
  lambda.space.input <- density(X, sigma=sigma.s)
  Y0 <- rpoint(n=n, f=lambda.space.input)
  
  # Initial configuration, temporal coordinates
  aux.prob <- pnorm((end.time-X$marks)/sigma.t) - pnorm((start.time-X$marks)/sigma.t)
  aux.components <- sample(x=1:n, size=n, replace=TRUE, prob=aux.prob)
  aux.marks <- rtruncnorm(n=n, mean=X$marks[aux.components], sd=sigma.t, a=start.time, b=end.time)
  Y <- Y0 %mark% aux.marks
  
  # Summaries of the initial configuration
  lambda.act <- intensityST.sep.grid(Y, sigma.t=sigma.t, sigma.s=sigma.s, times=ts, start.time=start.time, end.time=end.time)
  aux.lambda <- intensityST.sep.points(Y, sigma.t=sigma.t, sigma.s=sigma.s, start.time=start.time, end.time=end.time)
  L.act <- Lrt(Y, lambda=aux.lambda, Tmax=Tmax, Rmax=Rmax, start.time=start.time, end.time=end.time, grid=grid, gW=gW)
  
  # Check if L.act has NA values
  if (any(is.na(L.act))) {
    stop("NA values detected in L.act")
  }
  
  # Energy of the initial configuration
  en[,1] <- energy(lambda.target=lambda.input, L.target=L.input, lambda.act=lambda.act, L.act=L.act)[1]
  
  cat(paste0("Maximum number of iterations: ", max.it, ". Maximum number of rejected proposals in a row: ", max.reject, ".\n"))
  cat("Starting iterations, the ''+'' symbol indicates accepted proposal, the  ''.'' symbol indicates rejected proposal.\n\n")
  
  # Iteration step
  while (count < max.it & same <= max.reject){
    
    # New proposal
    i <- sample(1:n, 1)
    y <- rpoint(n=1, f=lambda.space.input)
    Ynew <- Y
    Ynew$x[i] <- y$x
    Ynew$y[i] <- y$y
    aux.comp <- sample(x=1:n, size=1, replace=TRUE, prob=aux.prob)
    aux.mark <- rtruncnorm(n=1, mean=X$marks[aux.comp], sd=sigma.t, a=start.time, b=end.time)
    Ynew$marks[i] <- aux.mark
    
    # Summaries of the proposed pattern
    lambda.act <- intensityST.sep.grid(Ynew, sigma.t=sigma.t, sigma.s=sigma.s, times=ts, start.time=start.time, end.time=end.time)
    aux.lambda <- intensityST.sep.points(Ynew, sigma.t=sigma.t, sigma.s=sigma.s, start.time=start.time, end.time=end.time)
    L.act <- Lrt(Ynew, lambda=aux.lambda, Tmax=Tmax, Rmax=Rmax, start.time=start.time, end.time=end.time, grid=grid, gW=gW)
    
    # Energy of the proposed pattern
    ennew <- energy(lambda.target=lambda.input, L.target=L.input, lambda.act=lambda.act, L.act=L.act)[1]
    
    # Accept or reject the proposal based on the energy
    if (ennew[1] <= en[1,count + 1]){
      Y <- Ynew
      en[,count + 2] <- ennew
      same <- 0} 
    else {
      en[,count + 2] <- en[,count + 1]
      same <- same + 1
    }
    count <- count + 1
    if (same==0){
      cat("+") # proposal accepted
    } else {
      cat(".") # proposal rejected
    }
    
    if ((count %% 100) == 0){
      cat("(finished ")
      cat(count)
      cat(" iterations so far)\n")
      par(mfrow=c(1,2))
      plot(en[1,1:(count+1)], type="l", main="Linear scale", ylab="Energy", xlab="Iterations", ylim=c(0, en[1,1]))
      points(x=(1:(count+1)), y=en[2,1:(count+1)], type="l", col=2) # L
      points(x=(1:(count+1)), y=en[3,1:(count+1)], type="l", col=3) # intensity
      legend("topright", legend=c("Total", "L", "Int"), lwd=c(1, 1, 1), col=c(1, 2, 3))
      plot(en[1,1:(count+1)], type="l", main="Logarithmic scale", log="y", ylab="Energy", xlab="Iterations", ylim=c(min(en[,count+1]), en[1,1]))
      points(x=(1:(count+1)), y=en[2,1:(count+1)], type="l", col=2) # L
      points(x=(1:(count+1)), y=en[3,1:(count+1)], type="l", col=3) # intensity
      legend("bottomleft", legend=c("Total", "L", "Int"), lwd=c(1, 1, 1), col=c(1, 2, 3))
    }
  }
  
  cat("\n")
  
  # Make sure the output pattern has the same window as the input pattern
  Y$window <- X$window
  
  # Determine the contribution of individual parts to the energy functional, useful for tuning the weights
  energy.final <- energy(lambda.target=lambda.input, L.target=L.input, lambda.act=lambda.act, L.act=L.act)
  EL.final <- energy.final[2]
  Elambda.final <- energy.final[3]
  
  # Create list of objects to be returned
  l = list(Output = Y, Count = count, Energy = en[,1:(count+1)], EL = EL.final, Elambda = Elambda.final)  
  return(l)
}
