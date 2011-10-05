plot.fptl <-
function (x, sfptl, from.t0 = TRUE, to.T = TRUE, dp.legend = TRUE, 
    dp.legend.cex = 1, ylab = TRUE, growth.points = TRUE, instants = TRUE, 
    ...) 
{
    if (!is.fptl(x)) 
        stop(paste(sQuote("x"), "is not of class", shQuote("fptl")))
    if (missing(sfptl)) {
        G <- growth.intervals(x$x, x$y)
        if (is.null(G)) {
            growth.points = FALSE
            instants = FALSE
            from.t0 = TRUE
            to.T = TRUE
        }
        else sfptl <- summary(x)
    }
    else {
        if (!is.summary.fptl(sfptl)) 
            stop(paste(sQuote("sfptl"), " object is not of class ", shQuote("summary.fptl")))
	
	  args <- as.list(attr(sfptl, "Call"))
	  if (is.element("zeroSlope", names(args))) zeroSlope <- eval(args$zeroSlope) else zeroSlope <- 0.01
	  if (is.element("p0.tol", names(args))) p0.tol <- eval(args$p0.tol) else p0.tol <- 8
	  if (is.element("k", names(args))) k <- eval(args$k) else k <- 3
        Y <- summary(x, zeroSlope, p0.tol, k)
        attr(Y, "Call") <- attr(sfptl, "Call")
        if (!identical(sfptl, Y)) 
            stop("the x and sfptl objects do not match up")
    }

    par(mar = c(3 + 2 * instants, 3 + ylab, 2, 1.5) + 0.1, ...)

    if (dp.legend) {	  
        args <- as.list(attr(x, "Call"))
	  if (is.element("env", names(args))){
		if (is.call(args$env)) args <- c(args[3:6], unlist(as.list(args$env)[-1]))
		else if (length(args$env) > 0) args <- c(args[3:6], unlist(args$env)) else args <- args[3:6]
	  }
	  else args <- args[3:6]
        dp.labels <- vector("expression", 3)
        dp.labels[[1]] <- substitute(paste("   Diffusion process:    ", 
            list(group("{", list(X(t), ~t ~ paste(" in  ", group("[", 
                list(t0, T), "]"))), "}"), ~~~P(X(t0) == x0) == 
                1)), args)
	  
	  logic <- unlist(lapply(args, function(x) if (is.call(x)) return(TRUE) else if (is.numeric(x)) any(x < 0, format(x) != format(x, scientific=FALSE)) else return(FALSE))) 		
	  if (any(logic)) args[logic] <- lapply(args[logic], function(x) as.call(parse(text=paste("(",deparse(x),")",sep="")))[[1]])				

	  dp.labels[[2]] <- substitute(paste(list(A[1](x, t) == m, 
            phantom(i)), ~~A[2](x, t) == v), list(m = eval(parse(text = paste("substitute(", 
            gsub("*", "%.%", attr(x, "dp")$mean, fixed=TRUE), ", args)", sep = ""))), v = eval(parse(text = paste("substitute(", 
            gsub("*", "%.%", attr(x, "dp")$var, fixed=TRUE), ", args)", sep = "")))))
        dp.labels[[3]] <- substitute("   Boundary:")
        dp.labels[[4]] <- substitute(S(t) == s, list(s = eval(parse(text = paste("substitute(", 
            gsub("*", "%.%", args$S, fixed=TRUE), ", args)", sep = "")))))
	  dp.h <- strheight(dp.labels, units = "inch") * sqrt(par("cex") * 
            dp.legend.cex)/par("pin")[2]
    }
    else dp.h <- numeric(4)
    m <- max(dp.h[c(2, 4)])
    if (growth.points) 
        gp.h <- strheight(expression(t[1]), units = "inch") * 
            sqrt(par("cex"))/par("pin")[2]
    else gp.h <- numeric(1)
    if (instants) 
        t0new.h <- strheight(expression(t[1]^{
            "*"
        }), units = "inch") * sqrt(par("cex"))/par("pin")[2]
    else t0new.h <- numeric(1)
    if (!from.t0) {
        if (x$x[1] < sfptl[1, 1]) {
            lg <- (x$x >= sfptl[1, 1])
            x$x <- x$x[lg]
            x$y <- x$y[lg]
        }
    }
    if (!to.T) {
        if (x$x[length(x$x)] > sfptl[nrow(sfptl), 5]) {
            j <- which(x$x >= sfptl[nrow(sfptl), 5])[1]
            j <- j + which.min(x$y[j:length(x$y)]) - 1
            x$x <- x$x[1:j]
            x$y <- x$y[1:j]
        }
    }
    plot(x$x, x$y, xlab = "", ylab = "", type = "l", las = 1, 
        ylim = c(-0.6 * t0new.h, 1 + 0.5 * gp.h + (dp.h[1] + 
            m)), axes = FALSE)
    box()
    axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), las = 1, mgp = c(3, 
        0.35 + 0.4/sqrt(par("cex")), 0), tcl = -0.35)
    axis(1, mgp = c(3, 0.35, 0), tcl = -0.35)
    title(main = "First-Passage-Time Location Function Plot", 
        line = 1)
    if (instants) 
        title(xlab = "t", line = 4)
    else title(xlab = "t", line = 1.6)
    if (ylab) 
        title(ylab = "FPTL(t)", line = 3)
    if (dp.legend) {
        x.labels <- par("usr")[1] + c(0, 0.25, 0.5, 0.75) * (par("usr")[2] - 
            par("usr")[1])
        y.labels <- c(par("usr")[4] - 0.5 * dp.h[1], par("usr")[4] - 
            dp.h[1] - 0.5 * m)
        text(x.labels[c(1, 3)], y.labels[1], dp.labels[c(1, 3)], 
            adj = 0, cex = dp.legend.cex)
        text(x.labels[c(2, 4)], y.labels[2], dp.labels[c(2, 4)], 
            adj = 0.5, cex = dp.legend.cex)
        abline(h = par("usr")[4] - dp.h[1] - m)
        segments(mean(par("usr")[1:2]), par("usr")[4], mean(par("usr")[1:2]), 
            par("usr")[4] - dp.h[1] - m)
    }
    if (any(growth.points, instants)) {
        if (growth.points) {
            segments(sfptl[, 1], par("usr")[3], sfptl[, 1], par("usr")[4] - 
                dp.h[1] - m, col = "darkgray", lwd = 1)
            text(sfptl[, 1] - 0.1 * par("cxy")[1], 1 + 0.5 * 
                (par("usr")[4] - dp.h[1] - m - 1), parse(text = paste("t[~ ", 
                1:nrow(sfptl), "]", sep = "")), adj = 1, col = "darkgray")
        }
        if (instants) {
            x.t <- matrix(sfptl[, c(2, 3, 5)], ncol = 3)
            y <- matrix(attr(sfptl, "FPTLValues")[, c(2, 3, 5)], 
                ncol = 3)
            points(x.t, y, type = "p", cex = min(1, 1/par("cex")))
            segments(x.t, par("usr")[3], x.t, y, lty = 8, lwd = 1)
            segments(par("usr")[1], y, x.t, y, lty = 8, lwd = 1)
            text(sfptl[, 2] - 0.4 * par("cxy")[1] * par("cex"), 
                par("usr")[3]/2, parse(text = paste("t[~ ", 1:nrow(x.t), 
                  "]^{~ ", shQuote("*"), "}", sep = "")))
            mtext(parse(text = paste("t[list(~ max,", 1:nrow(sfptl), 
                ")]^{~ ", shQuote("-"), "}", sep = "")), side = 1, 
                line = 1.6, at = sfptl[, 3], adj = 0, ...)
            mtext(parse(text = paste("t[list(~ max,", 1:nrow(sfptl), 
                ")]^{~ ", shQuote("+"), "}", sep = "")), side = 1, 
                line = 3, at = sfptl[, 5], adj = 0, ...)
        }
    }
}

