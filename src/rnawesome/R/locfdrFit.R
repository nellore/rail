## locfdrFit():
## arguments: for description of arguments, see locfdr() - exactly the same.
## return:  see locfdr() for many of the return elements, but this function has additional ones:
## --yt: heights of pink histogram bars that appear on the plots (i.e., heights of alt. density's histogram)
## --x: locations of pink histogram bars that appear on the plots (locations of alt. density's histogram)
## --mlest.lo and mlest.hi: if the function outputs a warning message saying "please re-run with mlest parameters = ...", said parameters are returned in these variables
## --needsfix: should the function advise you to re-run with different mlest parameters, needsfix will be 1, otherwise needsfix = 0.
## --nulldens: y-values of estimated null distribution density
## --fulldens: y-values of estimated full (mixture) density

locfdrFit <-
function (zz, bre = 120, df = 7, pct = 0, pct0 = 1/4, nulltype = 1, 
    type = 0, plot = 1, mult, mlests, main = " ", sw = 0, verbose=T) 
{
	require(locfdr)
	locmle = locfdr:::locmle
	loccov2 = locfdr:::loccov2
	loccov = locfdr:::loccov 
	mlest.lo <- mlest.hi <- yt <- x <- NULL
	needsfix <- 0
    call = match.call()
    if (length(bre) > 1) {
        lo <- min(bre)
        up <- max(bre)
        bre <- length(bre)
    }
    else {
        if (length(pct) > 1) {
            lo <- pct[1]
            up <- pct[2]
        }
        else {
            if (pct == 0) {
                lo <- min(zz)
                up <- max(zz)
            }
            if (pct < 0) {
                med = median(zz)
                ra = med + (1 - pct) * (range(zz) - med)
                lo = ra[1]
                up = ra[2]
            }
            if (pct > 0) {
                v <- quantile(zz, c(pct, 1 - pct))
                lo <- v[1]
                up <- v[2]
            }
        }
    }
    zzz <- pmax(pmin(zz, up), lo)
    breaks <- seq(lo, up, length = bre)
    zh <- hist(zzz, breaks = breaks, plot = F)
    x <- (breaks[-1] + breaks[-length(breaks)])/2
    yall <- y <- zh$counts
    K <- length(y)
    N <- length(zz)
    if (pct > 0) {
        y[1] <- min(y[1], 1)
        y[K] <- min(y[K], 1)
    }
    if (type == 0) {
        X <- cbind(1, ns(x, df = df))
        f <- glm(y ~ ns(x, df = df), poisson)$fit
    }
    else {
        X <- cbind(1, poly(x, df = df))
        f <- glm(y ~ poly(x, df = df), poisson)$fit
    }
    fulldens <- f
    l <- log(f)
    Fl <- cumsum(f)
    Fr <- cumsum(rev(f))
    D <- (y - f)/(f + 1)^0.5
    D <- sum(D[2:(K - 1)]^2)/(K - 2 - df)
    if (D > 1.5) 
        if(verbose){warning(paste("f(z) misfit = ", round(D, 1), ".  Rerun with increased df", 
            sep = ""))}
    if (nulltype == 3) {
        fp0 = matrix(NA, 6, 4)
        colnames(fp0) = c("delta", "sigleft", "p0", "sigright")
    }
    else {
        fp0 = matrix(NA, 6, 3)
        colnames(fp0) = c("delta", "sigma", "p0")
    }
    rownames(fp0) = c("thest", "theSD", "mlest", "mleSD", "cmest", 
        "cmeSD")
    fp0["thest", 1:2] = c(0, 1)
    fp0["theSD", 1:2] = 0
    imax <- seq(l)[l == max(l)][1]
    xmax <- x[imax]
    if (length(pct0) == 1) {
        pctup <- 1 - pct0
        pctlo <- pct0
    }
    else {
        pctlo <- pct0[1]
        pctup <- pct0[2]
    }
    lo0 <- quantile(zz, pctlo)
    hi0 <- quantile(zz, pctup)
    nx <- length(x)
    i0 <- (1:nx)[x > lo0 & x < hi0]
    x0 <- x[i0]
    y0 <- l[i0]
    if (nulltype == 3) {
        X00 <- cbind((x0 - xmax)^2, pmax(x0 - xmax, 0)^2)
    }
    else {
        X00 <- cbind(x0 - xmax, (x0 - xmax)^2)
    }
    lr <- lm(y0 ~ X00)
    co <- lr$coef
    if (nulltype == 3) {
        cmerror = I(is.na(co[3]) | is.na(co[2]))
        if (!cmerror) 
            cmerror = I(co[2] >= 0 | co[2] + co[3] >= 0)
    }
    else {
        cmerror = is.na(co[3])
        if (!cmerror) 
            cmerror = I(co[3] >= 0)
    }
    if (cmerror) {
        if (nulltype == 3) 
            stop("CM estimation failed.  Rerun with nulltype = 1 or 2.")
        else if (nulltype == 2) 
            stop("CM estimation failed.  Rerun with nulltype = 1.")
        else {
            X0 <- cbind(1, x - xmax, (x - xmax)^2)
            warning("CM estimation failed, middle of histogram non-normal")
        }
    }
    else {
        if (nulltype == 3) {
            X0 <- cbind(1, (x - xmax)^2, pmax(x - xmax, 0)^2)
            sigs <- 1/sqrt(-2 * (c(co[2], co[2] + co[3])))
            fp0["cmest", c(1, 2, 4)] <- c(xmax, sigs)
        }
        else {
            X0 <- cbind(1, x - xmax, (x - xmax)^2)
            xmaxx <- -co[2]/(2 * co[3]) + xmax
            sighat <- 1/sqrt(-2 * co[3])
            fp0["cmest", 1:2] <- c(xmaxx, sighat)
        }
        l0 <- as.vector(X0 %*% co)
        f0 <- exp(l0)
        p0 <- sum(f0)/sum(f)
        f0 <- f0/p0
        fp0["cmest", 3] <- p0
    }
    b = 4.3 * exp(-0.26 * log(N, 10))
    if (missing(mlests)) {
        med = median(zz)
        sc = diff(quantile(zz)[c(2, 4)])/(2 * qnorm(0.75))
        mlests = locmle(zz, xlim = c(med, b * sc))
        if (N > 5e+05) {
            if(verbose){warning("length(zz) > 500,000: For ML estimation, a wider interval than optimal was used.  To use the optimal interval, rerun with mlests = c(", 
                mlests[1], ", ", b * mlests[2], ").\n", sep = "")}
            mlest.lo = mlests[1]
            mlest.hi = b * mlests[2]
            needsfix = 1
            mlests = locmle(zz, xlim = c(med, sc))
        }
    }
    if (!is.na(mlests[1])) {
        if (N > 5e+05) 
            b = 1
        if (nulltype == 1) {
            Cov.in = list(x = x, X = X, f = f, sw = sw)
            ml.out = locmle(zz, xlim = c(mlests[1], b * mlests[2]), 
                d = mlests[1], s = mlests[2], Cov.in = Cov.in)
            mlests = ml.out$mle
        }
        else mlests = locmle(zz, xlim = c(mlests[1], b * mlests[2]), 
            d = mlests[1], s = mlests[2])
        fp0["mlest", 1:3] = mlests[1:3]
        fp0["mleSD", 1:3] = mlests[4:6]
    }
    if (sum(is.na(fp0[c(3, 5), 1:2])) == 0 & nulltype > 1) 
        if (abs(fp0["cmest", 1] - mlests[1]) > 0.05 | abs(log(fp0["cmest", 
            2]/mlests[2])) > 0.05) 
            warning("Discrepancy between central matching and maximum likelihood estimates.\nConsider rerunning with nulltype = 1")
    if (is.na(mlests[1])) {
        if (nulltype == 1) {
            if (is.na(fp0["cmest", 1])) 
                stop("CM and ML Estimation failed, middle of histogram non-normal")
            else stop("ML estimation failed.  Rerun with nulltype=2")
        }
        else warning("ML Estimation failed")
    }
    if (nulltype < 2) {
        delhat = xmax = xmaxx = mlests[1]
        sighat = mlests[2]
        p0 = mlests[3]
        f0 = dnorm(x, delhat, sighat)
        f0 = (sum(f) * f0)/sum(f0)
    }
    fdr = pmin((p0 * f0)/f, 1)
    f00 <- exp(-x^2/2)
    f00 <- (f00 * sum(f))/sum(f00)
    p0theo <- sum(f[i0])/sum(f00[i0])
    fp0["thest", 3] = p0theo
    fdr0 <- pmin((p0theo * f00)/f, 1)
    f0p <- p0 * f0
    if (nulltype == 0) 
        f0p <- p0theo * f00
    F0l <- cumsum(f0p)
    F0r <- cumsum(rev(f0p))
    Fdrl <- F0l/Fl
    Fdrr <- rev(F0r/Fr)
    Int <- (1 - fdr) * f * (fdr < 0.9)
    if (sum(x <= xmax & fdr == 1) > 0) 
        xxlo <- min(x[x <= xmax & fdr == 1])
    else xxlo = xmax
    if (sum(x >= xmax & fdr == 1) > 0) 
        xxhi <- max(x[x >= xmax & fdr == 1])
    else xxhi = xmax
    if (sum(x >= xxlo & x <= xxhi) > 0) 
        fdr[x >= xxlo & x <= xxhi] <- 1
    if (sum(x <= xmax & fdr0 == 1) > 0) 
        xxlo <- min(x[x <= xmax & fdr0 == 1])
    else xxlo = xmax
    if (sum(x >= xmax & fdr0 == 1) > 0) 
        xxhi <- max(x[x >= xmax & fdr0 == 1])
    else xxhi = xmax
    if (sum(x >= xxlo & x <= xxhi) > 0) 
        fdr0[x >= xxlo & x <= xxhi] <- 1
    if (nulltype == 1) {
        fdr[x >= mlests[1] - mlests[2] & x <= mlests[1] + mlests[2]] = 1
        fdr0[x >= mlests[1] - mlests[2] & x <= mlests[1] + mlests[2]] = 1
    }
    p1 <- sum((1 - fdr) * f)/N
    p1theo <- sum((1 - fdr0) * f)/N
    fall <- f + (yall - y)
    Efdr <- sum((1 - fdr) * fdr * fall)/sum((1 - fdr) * fall)
    Efdrtheo <- sum((1 - fdr0) * fdr0 * fall)/sum((1 - fdr0) * 
        fall)
    iup <- (1:K)[x >= xmax]
    ido <- (1:K)[x <= xmax]
    Eleft <- sum((1 - fdr[ido]) * fdr[ido] * fall[ido])/sum((1 - 
        fdr[ido]) * fall[ido])
    Eleft0 <- sum((1 - fdr0[ido]) * fdr0[ido] * fall[ido])/sum((1 - 
        fdr0[ido]) * fall[ido])
    Eright <- sum((1 - fdr[iup]) * fdr[iup] * fall[iup])/sum((1 - 
        fdr[iup]) * fall[iup])
    Eright0 <- sum((1 - fdr0[iup]) * fdr0[iup] * fall[iup])/sum((1 - 
        fdr0[iup]) * fall[iup])
    Efdr <- c(Efdr, Eleft, Eright, Efdrtheo, Eleft0, Eright0)
    Efdr[which(is.na(Efdr))] = 1
    names(Efdr) <- c("Efdr", "Eleft", "Eright", "Efdrtheo", "Eleft0", 
        "Eright0")
    if (nulltype == 0) 
        f1 <- (1 - fdr0) * fall
    else f1 <- (1 - fdr) * fall
    if (!missing(mult)) {
        mul = c(1, mult)
        EE = rep(0, length(mul))
        for (m in 1:length(EE)) {
            xe = sqrt(mul[m]) * x
            f1e = approx(xe, f1, x, rule = 2, ties = mean)$y
            f1e = (f1e * sum(f1))/sum(f1e)
            f0e = f0
            p0e = p0
            if (nulltype == 0) {
                f0e = f00
                p0e = p0theo
            }
            fdre = (p0e * f0e)/(p0e * f0e + f1e)
            EE[m] = sum(f1e * fdre)/sum(f1e)
        }
        EE = EE/EE[1]
        names(EE) = mul
    }
    Cov2.out = loccov2(X, X0, i0, f, fp0["cmest", ], N)
    Cov0.out = loccov2(X, matrix(1, length(x), 1), i0, f, fp0["thest", 
        ], N)
    if (sw == 3) {
        if (nulltype == 0) 
            Ilfdr = Cov0.out$Ilfdr
        else if (nulltype == 1) 
            Ilfdr = ml.out$Ilfdr
        else if (nulltype == 2) 
            Ilfdr = Cov2.out$Ilfdr
        else stop("With sw=3, nulltype must equal 0, 1, or 2.")
        return(Ilfdr)
    }
    if (nulltype == 0) 
        Cov = Cov0.out$Cov
    else if (nulltype == 1) 
        Cov = ml.out$Cov.lfdr
    else Cov = Cov2.out$Cov
    lfdrse <- diag(Cov)^0.5
    fp0["cmeSD", 1:3] = Cov2.out$stdev[c(2, 3, 1)]
    if (nulltype == 3) 
        fp0["cmeSD", 4] = fp0["cmeSD", 2]
    fp0["theSD", 3] = Cov0.out$stdev[1]
    if (sw == 2) {
        if (nulltype == 0) {
            pds = fp0["thest", c(3, 1, 2)]
            stdev = fp0["theSD", c(3, 1, 2)]
            pds. = t(Cov0.out$pds.)
        }
        else if (nulltype == 1) {
            pds = fp0["mlest", c(3, 1, 2)]
            stdev = fp0["mleSD", c(3, 1, 2)]
            pds. = t(ml.out$pds.)
        }
        else if (nulltype == 2) {
            pds = fp0["cmest", c(3, 1, 2)]
            stdev = fp0["cmeSD", c(3, 1, 2)]
            pds. = t(Cov2.out$pds.)
        }
        else stop("With sw=2, nulltype must equal 0, 1, or 2.")
        colnames(pds.) = names(pds) = c("p0", "delhat", "sighat")
        names(stdev) = c("sdp0", "sddelhat", "sdsighat")
        return(list(pds = pds, x = x, f = f, pds. = pds., stdev = stdev))
    }
    p1 <- seq(0.01, 0.99, 0.01)
    cdf1 <- rep(0, 99)
    fd <- fdr
    if (nulltype == 0) 
        fd <- fdr0
    for (i in 1:99) cdf1[i] <- sum(f1[fd <= p1[i]])
    cdf1 <- cbind(p1, cdf1/cdf1[99])
    mat <- cbind(x, fdr, Fdrl, Fdrr, f, f0, f00, fdr0, yall, 
        lfdrse, f1)
    namat <- c("x", "fdr", "Fdrleft", "Fdrright", "f", "f0", 
        "f0theo", "fdrtheo", "counts", "lfdrse", "p1f1")
    if (nulltype == 0) 
        namat[c(3, 4, 10)] <- c("Fdrltheo", "Fdrrtheo", "lfdrsetheo")
    dimnames(mat) <- list(NULL, namat)
    z.2 = rep(NA, 2)
    m = order(fd)[nx]
    if (fd[nx] < 0.2) {
        z.2[2] = approx(fd[m:nx], x[m:nx], 0.2, ties = mean)$y
    }
    if (fd[1] < 0.2) {
        z.2[1] = approx(fd[1:m], x[1:m], 0.2, ties = mean)$y
    }
    if(nulltype==0) nulldens <- p0theo*f00
    else nulldens <- p0*f0
    yt <- pmax(yall * (1 - fd), 0)
    if (plot > 0) {
        if (plot == 2 | plot == 3) 
            oldpar <- par(mfrow = c(1, 2), pty = "m")
        else if (plot == 4) 
            oldpar = par(mfrow = c(1, 3), pty = "m")
        hist(zzz, breaks = breaks, xlab = " ", main = main)
        yt <- pmax(yall * (1 - fd), 0)
        for (k in 1:K) lines(c(x[k], x[k]), c(0, yt[k]), lwd = 2, 
            col = 6)
        if (nulltype == 3) 
            title(xlab = paste("delta=", round(xmax, 3), "sigleft=", 
                round(sigs[1], 3), " sigright=", round(sigs[2], 
                  3), "p0=", round(fp0["cmest", 3], 3)))
        if (nulltype == 1 | nulltype == 2) 
            title(xlab = paste("MLE: delta:", round(mlests[1], 
                3), "sigma:", round(mlests[2], 3), "p0:", round(mlests[3], 
                3)), sub = paste("CME: delta:", round(fp0["cmest", 
                1], 3), "sigma:", round(fp0["cmest", 2], 3), 
                "p0:", round(fp0["cmest", 3], 3)))
        lines(x, f, lwd = 3, col = 3)
        if (nulltype == 0)
            lines(x, p0theo * f00, lwd = 2, lty = 2, col = 4)
        else lines(x, p0 * f0, lwd = 2, lty = 2, col = 4)
        if (!is.na(z.2[2])) 
            points(z.2[2], -0.5, pch = 24, col = "red", bg = "yellow")
        if (!is.na(z.2[1])) 
            points(z.2[1], -0.5, pch = 24, col = "red", bg = "yellow")
        if (nulltype == 1 | nulltype == 2) 
            Ef <- Efdr[1]
        else if (nulltype == 0) 
            Ef <- Efdr[4]
        if (plot == 2 | plot == 4) {
            if (nulltype == 0) 
                fdd <- fdr0
            else fdd = fdr
            matplot(x, cbind(fdd, Fdrl, Fdrr), type = "l", lwd = 3, 
                xlab = " ", ylim = c(0, 1.1), main = "fdr (solid); Fdr's (dashed)")
            title(xlab = paste("Efdr= ", round(Ef, 3)))
            abline(0, 0, lty = 3, col = 2)
            lines(c(0, 0), c(0, 1), lty = 3, col = 2)
        }
        if (plot == 3 | plot == 4) {
            if (sum(is.na(cdf1[, 2])) == nrow(cdf1)) 
                warning("cdf1 not available")
            else {
                plot(cdf1[, 1], cdf1[, 2], type = "l", lwd = 3, 
                  xlab = "fdr level", ylim = c(0, 1), ylab = "f1 proportion < fdr level", 
                  main = "f1 cdf of estimated fdr")
                title(sub = paste("Efdr= ", round(Ef, 3)))
                lines(c(0.2, 0.2), c(0, cdf1[20, 2]), col = 4, 
                  lty = 2)
                lines(c(0, 0.2), rep(cdf1[20, 2], 2), col = 4, 
                  lty = 2)
                text(0.05, cdf1[20, 2], round(cdf1[20, 2], 2))
                abline(0, 0, col = 2)
                lines(c(0, 0), c(0, 1), col = 2)
            }
        }
        if (plot > 1) 
            par(oldpar)
    }
    if (nulltype == 0) {
        ffdr <- approx(x, fdr0, zz, rule = 2, ties = "ordered")$y
    }
    else ffdr <- approx(x, fdr, zz, rule = 2, ties = "ordered")$y
    vl = list(fdr = ffdr, fp0 = fp0, Efdr = Efdr, cdf1 = cdf1, 
        mat = mat, z.2 = z.2)
    if (!missing(mult)) 
        vl$mult = EE
    vl$call = call
    vl$yt = yt
    vl$x = x
    vl$mlest.lo = mlest.lo
    vl$mlest.hi = mlest.hi
    vl$needsfix = needsfix
    vl$nulldens = nulldens
    vl$fulldens = fulldens
    vl
}
