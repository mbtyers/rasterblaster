#'Plots a kriging result in color!
#'@description Pass in a kriging result, and watch the world explode into technicolor, just like for Dorothy.
#'@param krigePlot The kriging object to use, returned from \link[geoR]{krige.conv}
#'@param type The type of plot to produce.  Specifying \code{type="predict"} will produce a map of kriging predictions, specifying \code{type="se"} will produce a map of kriging standard errors, and specifying \code{type="lower95"} or \code{type="upper95"} will produce a map of lower or upper 95 percent CIs of kriging predictions.  Defaults to \code{"predict"}.
#'@param z Alternately, a vector of quantities can be used, with the same length as the prediction locations.  Defaults to \code{NULL}.
#'@param ramp The color ramp to use for plotting.  Currently, allowed ramps are \code{"grey"}, \code{"red"}, \code{"rainbow"}, \code{"stoplight"}, and \code{"heat"}.  Defaults to \code{"rainbow"}.
#'@param invert Whether to invert the color ramp.  Try it both ways, see if it looks better.
#'@param alpha Color transparency to use for plotting.  This ranges from 0 (completely transparent) to 1 (true color).  Defaults to 0.8.
#'@param dark Color saturation to use for plotting.  This ranges from 0 (black) to 1 (true color).  Defaults to 0.85.
#'@param contour Whether to add contour lines.  Defaults to \code{TRUE}.
#'@param ... Additional plotting parameters
#'@author Matt Tyers
#'@examples
#'data(krigePlot,locations)
#'
#'rasterblaster(krigePlot)
#'rasterblaster(krigePlot,invert=T)
#'rasterblaster(krigePlot,alpha=.4)
#'rasterblaster(krigePlot,alpha=.4,dark=.7)
#'rasterblaster(krigePlot,ramp="grey")
#'@import geoR
#'@export
rasterblaster <- function(krigePlot,type="predict",z=NULL,ramp="rainbow",invert=F,alpha=.8,dark=.85,contour=T,xlab="X coord",ylab="Y coord",...) {
  locations <- eval(parse(text=krigePlot$call$locations))
  if(is.null(z)) {
    if(type=="predict") z <- krigePlot$predict
    if(type=="se") z <- sqrt(krigePlot$krige.var)
    if(type=="lower95") z <- krigePlot$predict-2*sqrt(krigePlot$krige.var)
    if(type=="upper95") z <- krigePlot$predict+2*sqrt(krigePlot$krige.var)
  }
  n <- dim(locations)[1]
  delx <- sort(unique(locations[,1]))[2] - sort(unique(locations[,1]))[1]
  dely <- sort(unique(locations[,2]))[2] - sort(unique(locations[,2]))[1]
  xleft <- locations[,1] - .5*delx
  xright <- locations[,1] + .5*delx
  ybottom <- locations[,2] - .5*dely
  ytop <- locations[,2] + .5*dely
  zcol <- (z-min(z))/(max(z)-min(z))
  if(!invert) zcol <- 1-zcol
  if(ramp=="grey") cols <- grey(zcol)
  if(ramp=="red") cols <- rgb(1,zcol,zcol)
  if(ramp=="rainbow") cols <- rainbow(n,start=0,end=.7)[rank(zcol)]
  if(ramp=="stoplight") cols <- rainbow(n,start=0,end=.3)[rank(zcol)]
  if(ramp=="heat") cols <- heat.colors(n)[rank(zcol)]
  cols <- adjustcolor(cols,alpha.f=alpha,red.f=dark,green.f=dark,blue.f=dark)
  plot(locations,asp=1,xlab=xlab,ylab=ylab,col="white",...=...)
  rect(xleft, ybottom, xright, ytop, col=cols, border=NA)
  if(contour) contour(krigePlot,val=z,add=T)
}

#'Makes a prediction grid from a geodata object
#'@description Shaboom!
#'@param geodata The geodata object to use
#'@param length The number of points to use in the X- and Y- direction.  Defaults to 100, resulting in a 100 x 100 grid.
#'@param resolution Alternately, the grid can be defined by spacing with this dimension.  Specifying \code{resolution=10} would result in a grid in which points are spaced 10 units apart in the X- and Y- directions.  Defaults to \code{NULL}.
#'@author Matt Tyers
#'@examples
#'data(lnLa)
#'
#'grid <- makegrid(lnLa)
#'plot(lnLa$coords)
#'points(grid,pch='.')
#'@import geoR
#'@export
makegrid <- function(geodata, length=100,resolution=NULL) {
  if(is.null(resolution)) {
    xx <- seq(min(geodata$coords[,1]),max(geodata$coords[,1]),length=length)
    yy <- seq(min(geodata$coords[,2]),max(geodata$coords[,2]),length=length)
  }
  if(!is.null(resolution)) {
    xx <- seq(min(geodata$coords[,1]),max(geodata$coords[,1]),by=resolution)
    yy <- seq(min(geodata$coords[,2]),max(geodata$coords[,2]),by=resolution)
  }
  gridd <- expand.grid(xx,yy)
  return(gridd)
}

#'Likelihood ratio test for two likfit models
#'@description Never have to remember which way to subtract again!
#'@param likfit1 The first model returned from \link[geoR]{likfit}
#'@param likfit1 The second model returned from \link[geoR]{likfit}.  It doesn't matter which is which.
#'@return The p-value from the associated likelihood ratio test
#'@note Don't even think about using non-nested models, or the wrath of Dr. Barry will haunt you for the rest of your days.
#'@author Matt Tyers
#'@examples
#'data(lnLa)
#'model1 <- likfit(lnLa,cov.model="exponential",ini.cov.pars=c(1,18))
#'model2 <- likfit(lnLa,cov.model="exponential",ini.cov.pars=c(1,18),trend="1st")
#'
#'logliktest(model1,model2)
#'logliktest(model2,model1)
#'@import geoR
#'@export
logliktest <- function(likfit1,likfit2) {
  if(likfit1$npars < likfit2$npars) pval <- 1-pchisq(2*(likfit2$loglik-likfit1$loglik), likfit2$npars-likfit1$npars)
  if(likfit1$npars > likfit2$npars) pval <- 1-pchisq(2*(likfit1$loglik-likfit2$loglik), likfit1$npars-likfit2$npars)
  return(pval)
}
