## getParams()
## arguments: 
## --tstats: vector of t statistics (probably the output of getTstats()) obtained by testing for DE at each base
## --plots: if TRUE, diagnostic plots are generated
## --plotfile: (optional) string containing a file where plots will be written, if plots = TRUE. If plots=TRUE and no plotfile is provided, interactive graphics will be used.
## --verbose: if TRUE, messages will print as the function runs and parameters are estimated.
## return:
## a list containing several elements:
## $params: a list containing $mean and $sd, vectors giving the estimated mean and sd (respectively) of the zero state, null distribution, DE-up distribution, and DE-down distribution (respectively)
## $stateprobs: a vector containing the estimated proportion of zeros, null bases, DE-up bases, and DE-down bases, respectively
## $pi.z, $pi.0, and $pi.1: just provides easier access to stateprobs: note that stateprobs = c(pi.z, pi.0, pi.1/2, pi.1/2)
## suggested use: findparams.chr1 = getParams(tstats); regions = getRegions("HMM",1,pos,tstats,stateprobs=findparams.chr1$stateprobs, params=findparams.chr1$params)
## (i.e., params and stateprobs elements correspond to what should be fed to getRegions()


getParams = function(tstats, plots = FALSE, plotfile = NULL, verbose = F){
          if (plots & !is.null(plotfile)) pdf(file = plotfile)
          pi.z = sum(tstats == 0)/length(tstats)
          nzt = tstats[tstats != 0]
          if (plots) {
              hist(nzt, freq = FALSE, breaks = 100, xlab = "t", main = "Histogram of nonzero t-statistics")
              zzaxis = seq(-5, 5, length = 1000)
              lines(zzaxis, dnorm(zzaxis), col = "red")
              legend("topright", "N(0,1)", lty = 1, col = "red", bty = "n")
          }
          plotarg = ifelse(plots, 1, 0)
          fdrmodel.defaults = try(locfdrFit(nzt, plot = plotarg, verbose = F, main = "initial locfdr",nulltype=0))
          if(class(fdrmodel.defaults)=="try-error"){
               params = list(mean = c(0,1,2,-2), sd = c(0.0000001,1,1,1))
               pi.0 = 0.99*(1-pi.z)
               pi.1 = 1-pi.z-pi.0
               out = list(params = params, stateprobs = c(pi.z, pi.0, pi.1/2, pi.1/2))
               warning("locfdr estimation failed. Defaulting to N(0,1) null and N(+/-2, 1) alternatives")
               return(out)
          }
          if (fdrmodel.defaults$needsfix == 1) {
             fdrmodel <- locfdrFit(nzt, mlest = c(fdrmodel.defaults$mlest.lo,fdrmodel.defaults$mlest.hi), verbose = F, plot = plotarg,
                   main = "locfdr, incorporating mlest",nulltype=0)
             if (fdrmodel$needsfix == 1) stop("problem with mlest parameters in locfdr")
          }
          if (fdrmodel.defaults$needsfix == 0) fdrmodel <- fdrmodel.defaults
          pi.0.tmp = fdrmodel$fp0[3, 3]
          if (pi.0.tmp >= 1){
               warning("Proportion of null nonzero bases estimated to be greater than 1.  Defaulting to 0.999.")
               pi.0.tmp <- 0.999
          }
          pi.0 = pi.0.tmp * (1 - pi.z)
          pi.1 = 1 - pi.z - pi.0
          if (verbose) message("starting DEup.mean...")
          u = find.mean.up(0.99, null.mean = fdrmodel$fp0[3, 1], null.sd = fdrmodel$fp0[3,2], null.prop = pi.0.tmp, vals = nzt)
          DEup.mean = u$m
          if (verbose) message("starting DEdown.mean...")
          d = find.mean.down(0.01, null.mean = fdrmodel$fp0[3, 1],null.sd = fdrmodel$fp0[3, 2], null.prop = pi.0.tmp, vals = nzt)
          DEdown.mean = d$m
          if (names(u)[2] == "p") {
              DEup.sd = find.sd(prev.p = u$p, found.mean = DEup.mean, null.mean = fdrmodel$fp0[3, 1], null.sd = fdrmodel$fp0[3, 2], null.prop = pi.0.tmp, vals = nzt, up = TRUE)
          }
          if (names(u)[2] == "s") {
              DEup.sd = u$s
          }
          if (names(d)[2] == "p") {
              DEdown.sd = find.sd(prev.p = d$p, found.mean = DEdown.mean,
                       null.mean = fdrmodel$fp0[3, 1], null.sd = fdrmodel$fp0[3, 2], null.prop = pi.0.tmp, vals = nzt, up = FALSE)
          }
          if (names(d)[2] == "s") {
                DEdown.sd = d$s
          }
         if (plots) {
              zzaxis <- seq(-10, 10, length = 2000)
              hist(tstats, breaks = 100, col = "gray", freq = F, ylim = c(0, 0.4), main = "fitted distributions", xlim =  c(-5, 5), xlab = "t statistics")
              lines(zzaxis, pi.0 * dnorm(zzaxis, fdrmodel$fp0[3, 1],fdrmodel$fp0[3, 2]), lwd = 3)
              lines(zzaxis, (pi.1)/2 * dnorm(zzaxis, DEup.mean, DEup.sd),lwd = 3, col = "red")
              lines(zzaxis, (pi.1)/2 * dnorm(zzaxis, DEdown.mean, DEdown.sd), lwd = 3, col = "green")
              legend("topright", lty = c(1, 1, 1), lwd = c(3, 3, 3), col = c("black", "red", "green"), c("null", "DE up", "DE down"))
              if (!is.null(plotfile))  dev.off()
          }
          params = list(mean = c(0, fdrmodel$fp0[3, 1], DEup.mean, DEdown.mean), sd = c(1e-07, fdrmodel$fp0[3, 2], DEup.sd, DEdown.sd))
          out = list(params = params, stateprobs = c(pi.z, pi.0, pi.1/2, pi.1/2))
          return(out)
        }
         
