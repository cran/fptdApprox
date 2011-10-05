plot.fpt.density <-
function (x, from.t0, to.T, dp.legend = TRUE, dp.legend.cex = 1, 
    ylab = TRUE, growth.points = FALSE, instants = FALSE, ...) 
{
    if (!is.fpt.density(x)) 
        stop(paste(sQuote("x"), "is not of class", shQuote("fpt.density")))

    args <- as.list(attr(x, "Call"))
    if (is.element("from.t0", names(args))) fromt0Call <- args$from.t0 else fromt0Call <- FALSE

    if (missing(from.t0)) 
        from.t0 <- fromt0Call

    if (is.element("to.T", names(args))) toTCall <- args$to.T else toTCall <- FALSE

    if (missing(to.T)) 
        to.T <- toTCall

    args <- as.list(attr(attr(x, "summary.fptl"), "FPTLCall"))
    if (is.element("env", names(args))){
	if (is.call(args$env)) args <- c(args[3:6], unlist(as.list(args$env)[-1]))
	else if (length(args$env) > 0) args <- c(args[3:6], unlist(args$env)) else args <- args[3:6]
    }
    else args <- args[3:6]
    
    par(mar = c(3 + 2 * instants, 3 + ylab, 2, 1.5) + 0.1, ...)

    if (dp.legend) {	          
	  logic <- unlist(lapply(args, function(x) if (is.call(x)) return(TRUE) else if (is.numeric(x)) any(x < 0, format(x) != format(x, scientific=FALSE)) else return(FALSE))) 		
	  if (any(logic)) args[logic] <- lapply(args[logic], function(x) as.call(parse(text=paste("(",deparse(x),")",sep="")))[[1]])		
		
        dp.labels <- vector("expression", 3)
        dp.labels[[1]] <- substitute(paste("   Diffusion process:    ", 
            list(group("{", list(X(t), ~t ~ paste(" in  ", group("[", 
                list(t0, T), "]"))), "}"), ~~~P(X(t0) == x0) == 
                1)), args)
	  dp.labels[[2]] <- substitute(paste(list(A[1](x, t) == m, 
            phantom(i)), ~~A[2](x, t) == v), list(m = eval(parse(text = paste("substitute(", 
            gsub("*", "%.%", attr(attr(x, "summary.fptl"), "dp")$mean, fixed=TRUE), ", args)", sep = ""))), v = eval(parse(text = paste("substitute(", 
            gsub("*", "%.%", attr(attr(x, "summary.fptl"), "dp")$var, fixed=TRUE), ", args)", sep = "")))))
        dp.labels[[3]] <- substitute("   Boundary:")
        dp.labels[[4]] <- substitute(S(t) == s, list(s = eval(parse(text = paste("substitute(", 
            gsub("*", "%.%", args$S, fixed=TRUE), ", args)", sep = "")))))
	  dp.h <- strheight(dp.labels, units = "inch") * sqrt(par("cex") * 
            dp.legend.cex)/par("pin")[2]
    }
    else dp.h <- numeric(4)
    
    args <- lapply(args, eval)

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
    if ((!from.t0) & fromt0Call & (x$x[1] < attr(x, "summary.fptl")[1, 
        1])) {
        lg <- (x$x >= attr(x, "summary.fptl")[1, 1])
        x$x <- x$x[lg]
        x$y <- x$y[lg]
    }
    if (from.t0 & (!fromt0Call) & (args$t0 < x$x[1])) {
        x$x <- c(args$t0, x$x)
        x$y <- c(0, x$y)
    }
    if ((!to.T) & toTCall & (attr(x, "Steps")[nrow(attr(x, "Steps")), 
        2] > attr(x, "summary.fptl")[length(attr(x, "summary.fptl"))])) {
        lg <- (x$x <= attr(x, "summary.fptl")[length(attr(x, 
            "summary.fptl"))])
        x$x <- x$x[lg]
        x$y <- x$y[lg]
    }
    ymax <- max(x$y)
    plot(x$x, x$y, xlab = "", ylab = "", type = "l", las = 1, 
        ylim = c(-0.6 * t0new.h * ymax, ymax + (0.5 * gp.h + 
            dp.h[1] + m) * ymax), axes = FALSE)
    box()
    ticks <- axTicks(2)
    axis(2, las = 1, at = ticks[ticks <= par("usr")[4] - (dp.h[1] + 
        m) * ymax], mgp = c(3, 0.35 + 0.4/sqrt(par("cex")), 0), 
        tcl = -0.35)
    axis(1, mgp = c(3, 0.35, 0), tcl = -0.35)
    title(main = "Approximate First-Passage-Time Density Function Plot", 
        line = 1)
    if (instants) 
        title(xlab = "t", line = 4)
    else title(xlab = "t", line = 1.6)
    if (ylab) 
        title(ylab = parse(text = "g[1](t)"), line = 3)
    if (dp.legend) {
        x.labels <- par("usr")[1] + c(0, 0.25, 0.5, 0.75) * (par("usr")[2] - 
            par("usr")[1])
        y.labels <- c(par("usr")[4] - 0.5 * dp.h[1] * ymax, par("usr")[4] - 
            (dp.h[1] + 0.5 * m) * ymax)
        text(x.labels[c(1, 3)], y.labels[1], dp.labels[c(1, 3)], 
            adj = 0, cex = dp.legend.cex)
        text(x.labels[c(2, 4)], y.labels[2], dp.labels[c(2, 4)], 
            adj = 0.5, cex = dp.legend.cex)
        abline(h = par("usr")[4] - (dp.h[1] + m) * ymax)
        segments(mean(par("usr")[1:2]), par("usr")[4], mean(par("usr")[1:2]), 
            par("usr")[4] - (dp.h[1] + m) * ymax)
    }
    if (any(growth.points, instants)) {
        Y <- attr(x, "summary.fptl")
        if (growth.points) {
            i <- which(Y[, 1] >= x$x[1])
            if (length(i) > 0) {
                segments(Y[i, 1], par("usr")[3], Y[i, 1], par("usr")[4] - 
                  (dp.h[1] + m) * ymax, col = "darkgray", lwd = 1)
                text(Y[i, 1] - 0.1 * par("cxy")[1], par("usr")[4] - 
                  (dp.h[1] + m) * ymax - 0.5 * gp.h * ymax, parse(text = paste("t[~ ", 
                  i, "]", sep = "")), adj = 1, col = "darkgray")
            }
        }
        if (instants) {
            x.t <- matrix(Y[, c(2, 3, 5)], ncol = 3)
            segments(x.t, par("usr")[3], x.t, par("usr")[4] - 
                (dp.h[1] + m) * ymax, lty = 8, lwd = 1)
            text(Y[, 2] - 0.4 * par("cxy")[1] * par("cex"), par("usr")[3]/2, 
                parse(text = paste("t[~ ", 1:nrow(x.t), "]^{~ ", 
                  shQuote("*"), "}", sep = "")))
            mtext(parse(text = paste("t[list(~ max,", 1:nrow(Y), 
                ")]^{~ ", shQuote("-"), "}", sep = "")), side = 1, 
                line = 1.6, at = Y[, 3], adj = 0, ...)
            mtext(parse(text = paste("t[list(~ max,", 1:nrow(Y), 
                ")]^{~ ", shQuote("+"), "}", sep = "")), side = 1, 
                line = 3, at = Y[, 5], adj = 0, ...)
        }
    }
}

