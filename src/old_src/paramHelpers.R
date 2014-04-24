# helper functions for getParams.R
# 7/18/12
# updated 8/7/12 to include the failsafe function
# updated again 8/27/12 - new failsafe.

# little helper function
last <- function(x) return(tail(x, n=1))

## getParams.failsafe:
## arguments: null.mean and null.sd, mean and sd of estimated null distribution
## return: list with $DEup.mean, $DEdown.mean, $DEup.sd, and $DEdown.sd, giving parameters of alternative distributions as estimated by locfdr
## used in case our numerical estimation method fails.
## just use the 5th (for DE down) and 95th (for DE up) percentiles of the null distribution as alternative means
## and then have the alternative SDs be equal to the null SD.
getParams.failsafe <- function(null.mean, null.sd){
		DEup.mean = qnorm(.95, mean=null.mean, sd=null.sd)
		DEup.sd = DEdown.sd = null.sd
		DEdown.mean = qnorm(.05, mean=null.mean, sd=null.sd)
		return(list(DEup.mean = DEup.mean, DEdown.mean = DEdown.mean, DEup.sd = DEup.sd, DEdown.sd = DEdown.sd))
}

## get.numalts():
## arguments:
## --pctil: a percentile (e.g. 0.01 for first percentile, 0.88 for 88th percentile...)
## --null.mean, null.sd: parameters of null distribution as estimated from locfdr
## --null.prop: proportion of null values as estimated from locfdr
## --vals: full vector of t statistics (or whatever) that you estimated your null/alt distributions from
## --up: TRUE if you're estimating the DEup parameters, FALSE if estimating DEdown parameters.
## return: the estimated number of alternative values (in vals) above (below, if up=F) the pctil'th percentile of the null distribution
## the above number is in element $num, the element $val contains the pctil'th percentile of the null distribution.
## note that this function probably shouldn't be used directly: only used by find.mean and find.sd (below)
get.numalts <- function(pctil,null.mean,null.sd,null.prop,vals,up = TRUE){
	cutoff.val = qnorm(pctil,null.mean,null.sd)
	num.nulls = null.prop*length(vals)
	# I use "above" to mean "above or below". very loose...
	coef1 = ifelse(up,1-pctil,pctil)
	num.nulls.above = coef1*num.nulls
	num.total.above = ifelse(up,sum(vals>cutoff.val),sum(vals<cutoff.val))
	num.alts.above = num.total.above - num.nulls.above
	if(num.alts.above<=0) stop("negative num.alts.above. rethink algorithm. sad day.")
	return(list(num=round(num.alts.above),val=cutoff.val))
}

## find.mean.up():
## arguments:
## --init.value: intial percentile of null distribution we are going to look at
## --null.mean, null.sd: parameters of null distribution as estimated from locfdr
## --null.prop: proportion of null values as estimated from locfdr
## --vals: full vector of t statistics (or whatever) that you estimated your null/alt distributions from
## return:
## --if our method SUCCEEDS: list with elements $m (the estimated mean of the DE-up distribution) and $p (the percentile of the null distribution used to calculate that mean)
## --if our method FAILS: list with elements $m (estimated mean of the DE-up distribution) and $s (estimated sd of DE-up dist), as estimated with getParams.failsafe().
find.mean.up <- function(init.value, null.mean, null.sd, null.prop, vals){
	if(init.value<=0.5) stop("finding DE up quantile - init.value should be >0.5")
	x = init.value
	target = round(0.5*0.5*(1-null.prop)*length(vals))  #half the expected # of DEup values.
	history = pseq = numseq = qseq = NULL
	counter = 0
	while(TRUE){
		iter = try(get.numalts(x,null.mean,null.sd,null.prop,vals,up=T),silent=T)
		if(class(iter)=="try-error") {
			warning("Numerical estimation of DE-up mean failed.  Defaulting to mean = 95th percentile of estimated null distribution, sd = sd of estimated null distribution.")
			gp = getParams.failsafe(null.mean, null.sd)
			return(list(m=gp$DEup.mean, s=gp$DEup.sd))
		}
		if(iter$num == target) return(list(m=iter$val,p = pseq[counter]))
		if(counter>1){
			if(abs(pseq[counter]-pseq[counter-1])<1e-8) return(list(m=iter$val,p=pseq[counter]))
		}
		counter = counter+1
		numseq[counter] = iter$num
		qseq[counter] = iter$val
		pseq[counter] = x

		if(counter>1){
			if(history[counter-1]==1 & numseq[counter]>=numseq[counter-1]){
				warning("Numerical estimation of DE-up mean failed.  Defaulting to mean = 95th percentile of estimated null distribution, sd = sd of estimated null distribution.")
				gp = getParams.failsafe(null.mean, null.sd)
				return(list(m=gp$DEup.mean, s=gp$DEup.sd))
			}			
			if(history[counter-1]==-1 & numseq[counter]<=numseq[counter-1]){
				warning("Numerical estimation of DE-up mean failed.  Defaulting to mean = 95th percentile of estimated null distribution, sd = sd of estimated null distribution.")
				gp = getParams.failsafe(null.mean, null.sd)
				return(list(m=gp$DEup.mean, s=gp$DEup.sd))
			}
		}

		if(iter$num > target){
			history[counter] = 1
			if(counter==1){ x = (x+1)/2; next }
		}
		if(iter$num < target){
			history[counter] = -1
			if(counter==1){ x = (x+0.5)/2; next }
		}

		if(history[counter-1] != history[counter]) x = (x+pseq[counter-1])/2
		if(history[counter-1] == history[counter]){
				prev = last(which(history != history[counter]))
				if(length(prev)==0 & iter$num > target){x = (x+1)/2; next}
				if(length(prev)==0 & iter$num < target){x = (x+0.5)/2; next}
				x = (x+pseq[prev])/2
			}
		}
}


## find.mean.down():  slightly modified version of find.mean.up...
## arguments:
## --init.value: intial percentile of null distribution we are going to look at
## --null.mean, null.sd: parameters of null distribution as estimated from locfdr
## --null.prop: proportion of null values as estimated from locfdr
## --vals: full vector of t statistics (or whatever) that you estimated your null/alt distributions from
## return:
## --if our method SUCCEEDS: list with elements $m (the estimated mean of the DE-down distribution) and $p (the percentile of the null distribution used to calculate that mean)
## --if our method FAILS: list with elements $m (estimated mean of the DE-down distribution) and $s (estimated sd of DE-down dist), as estimated with getParams.failsafe().
find.mean.down <- function(init.value, null.mean, null.sd, null.prop, vals){
	if(init.value>=0.5) stop("finding DE down quantile - init.value should be <0.5")
	x = init.value
	target = round(0.5*0.5*(1-null.prop)*length(vals))  #half the expected # of DEdown values.
	history = pseq = numseq = qseq = NULL
	counter = 0
	while(TRUE){
		iter = try(get.numalts(x,null.mean,null.sd,null.prop,vals,up=F),silent=T)
		if(class(iter)=="try-error") {
			warning("Numerical estimation of DE-down mean failed.  Defaulting to mean = 5th percentile of estimated null distribution, sd = sd of estimated null distribution.")
			gp = getParams.failsafe(null.mean, null.sd)
			return(list(m=gp$DEdown.mean,s=gp$DEdown.sd))
		}
		if(iter$num == target) return(list(m=iter$val,p = pseq[counter]))
		if(counter>1){
			if(abs(pseq[counter]-pseq[counter-1])<1e-8) return(list(m=iter$val,p=pseq[counter]))
		}
		counter = counter+1
		numseq[counter] = iter$num
		qseq[counter] = iter$val
		pseq[counter] = x

		if(counter>1){
			if(history[counter-1]==1 & numseq[counter]>=numseq[counter-1]){
				warning("Numerical estimation of DE-down mean failed.  Defaulting to mean = 5th percentile of estimated null distribution, sd = sd of estimated null distribution.")
				gp = getParams.failsafe(null.mean, null.sd)
				return(list(m=gp$DEdown.mean, s=gp$DEdown.sd))
			}
			if(history[counter-1]==-1 & numseq[counter]<=numseq[counter-1]){
				warning("Numerical estimation of DE-down mean failed.  Defaulting to mean = 5th percentile of estimated null distribution, sd = sd of estimated null distribution.")
				gp = getParams.failsafe(null.mean, null.sd)
				return(list(m=gp$DEdown.mean, s=gp$DEdown.sd))
			}
		}

		if(iter$num > target){
			history[counter] = 1
			if(counter==1){ x = x/2; next }
		}
		if(iter$num < target){
			history[counter] = -1
			if(counter==1){ x = (x+0.5)/2; next }
		}

		if(history[counter-1] != history[counter]) x = (x+pseq[counter-1])/2
		if(history[counter-1] == history[counter]){
				prev = last(which(history != history[counter]))
				if(length(prev)==0 & iter$num > target){x = x/2; next}
				if(length(prev)==0 & iter$num < target){x = (x+0.5)/2; next}
				x = (x+pseq[prev])/2
			}
		}
}

## just used in case you don't want to use the above find.mean.up or find.mean.down.  arguments/return are the same.
find.mean <- function(init.value, null.mean, null.sd, null.prop, vals, up = TRUE){
	if(up) k = find.mean.up(init.value, null.mean, null.sd, null.prop, vals)
	if(!up) k = find.mean.down(init.value, null.mean, null.sd, null.prop, vals)
	return(k)
}

## find.sd(): used to find SD of alternative distribution if find.mean.up or down succeeds
## arguments:
## --prev.p: percentile used to find the alternative mean (i.e., the $p return from find.mean.up or down)
## --found.mean: alternative mean found (i.e., the $m return from find.mean.up or down)
## --null.mean, null.sd: parameters of null distribution as estimated from locfdr
## --null.prop: proportion of null values as estimated from locfdr
## --vals: full vector of t statistics (or whatever) that you estimated your null/alt distributions from
## --up: TRUE if you want the DEup sd, FALSE if you want the DEdown.sd
## return:
## --if this method succeeds, the estimated standard deviation of the alternative distribution (using this method)
## --if this method succeeds, the estimated sd of the alternative distribution using the output from locfdr. warning message is printed.
find.sd <- function(prev.p, found.mean, null.mean, null.sd, null.prop, vals, up=T){
	if(up) new.percentile = (1+prev.p)/2 # we increase the percentile we're going to look at
	if(!up) new.percentile = prev.p/2 #decrease it.
	cutoff.val = qnorm(new.percentile, null.mean, null.sd) # x-value that is that percentile of the theoretical null distribution
	num.nulls = null.prop*length(vals)
	coef1 = ifelse(up,1-new.percentile,new.percentile)
	num.nulls.above = coef1*num.nulls
	num.total.above = ifelse(up,sum(vals>cutoff.val),sum(vals<cutoff.val))
	num.alts.above = num.total.above - num.nulls.above
	num.alts = (1-null.prop)*0.5*length(vals) # this is the number of DEup values
	if(up) alt.percentile = 1-(num.alts.above/num.alts)
	if(!up) alt.percentile = num.alts.above/num.alts
	zstat = qnorm(alt.percentile)
	if(is.na(zstat)){
		warning("Numerical standard deviation estimation failed.  Defaulting to sd of estimated null distribution.")
		gp = getParams.failsafe(null.mean, null.sd)
		if(up) return(gp$DEup.sd)
		if(!up) return(gp$DEdown.sd)
	}
	sigma = (cutoff.val - found.mean)/zstat
	if(num.alts.above<=0 | sigma<=0) {
		warning("Numerical standard deviation estimation failed. Defaulting to sd of estimated null distribution.")
		gp = getParams.failsafe(null.mean, null.sd)
		if(up) return(gp$DEup.sd)
		if(!up) return(gp$DEdown.sd)
	}
	return(sigma)
}
