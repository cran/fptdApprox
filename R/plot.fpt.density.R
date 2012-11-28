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
   
    par(cex = 1, ps = 9)      
    par(mar = c(3 + 2 * instants, 3.25 + ylab, 2.5, 1) + 0.1, ...)        
    pin.height <- par("pin")[2] - sum(par("mai")[c(1,3)]*(par("cex")-1))

    if (dp.legend) {
	  logic <- unlist(lapply(args, function(x) if (is.call(x)) return(TRUE) else if (is.numeric(x)) any(x < 0, format(x) != format(x, scientific=FALSE)) else return(FALSE))) 		
	  if (any(logic)) args[logic] <- lapply(args[logic], function(x) as.call(parse(text=paste("(",deparse(x),")",sep="")))[[1]])		
			          
        dp.labels <- vector("expression", 3)
        dp.labels[[1]] <- substitute(paste(" Diffusion process: ", 
            list(group("{", list(X(t), paste(t, "  in  ", group("[", 
                list(t0, T), "]"))), "}"), ~~P(X(t0) == x0) == 
                1)), args)
	  
	  dp.labels[[2]] <- substitute(list(paste(A[1](x, t) == m, " "), 
            ~~A[2](x, t) == v), list(m = eval(parse(text = paste("substitute(", 
            gsub("*", "%.%", attr(attr(x, "summary.fptl"), "dp")$mean, fixed=TRUE), ", args)", sep = ""))), v = eval(parse(text = paste("substitute(", 
            gsub("*", "%.%", attr(attr(x, "summary.fptl"), "dp")$var, fixed=TRUE), ", args)", sep = "")))))

        dp.labels[[3]] <- substitute(paste(" Boundary: ", S(t) == s), list(s = eval(parse(text = paste("substitute(", 
            gsub("*", "%.%", args$S, fixed=TRUE), ", args)", sep = "")))))
        dp.width <- sum(strwidth(dp.labels[1:2], units = "inch", cex = dp.legend.cex)) + strwidth(" ,  ", units = "inch", cex = dp.legend.cex)
	  pin.width <- par("pin")[1] - sum(par("mai")[c(2,4)]*(par("cex")-1))	  
	  dp.h <- strheight(dp.labels, units = "inch", cex = dp.legend.cex)
	  logic <- (dp.width < pin.width)
	  if (logic) dp.height <- max(dp.h[1:2]) + dp.h[3] + 0.3*par("cin")[2] else dp.height <- sum(dp.h) + 0.4*par("cin")[2]
    }
    else dp.height <- numeric(1)
    
    gp.h <- 0.03 * pin.height
    ti.h <- gp.h
    if (growth.points) gp.h <- max(gp.h, 1.25 * strheight(expression(t), units = "inch", cex = sqrt(par("cex"))) + 0.2*par("cin")[2])
    if (instants) ti.h <- max(ti.h, 1.5 * strheight(expression(t[1]), units = "inch", cex = sqrt(par("cex"))) + 0.2*par("cin")[2])
    
    args <- lapply(args, eval)

    if ((!from.t0) & fromt0Call & (x$x[1] < attr(x, "summary.fptl")[1,1])) {
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
        lg <- (x$x <= attr(x, "summary.fptl")[length(attr(x, "summary.fptl"))])
        x$x <- x$x[lg]
        x$y <- x$y[lg]
    }
    ymax <- max(x$y)

    h <- pin.height - dp.height - gp.h - ti.h
    A <- - ti.h * ymax/h
    B <- (1 + (dp.height + gp.h)/h) * ymax
    y2 <- ((B-A) + 1.08*(A+B))/2.16
    y1 <- A + B - y2
    plot(x$x, x$y, xlab = "", ylab = "", type = "l", las = 1, ylim = c(y1, y2), axes = FALSE)
    box()    

    y3 <- par("usr")[4]
    if (dp.legend) {
	  dp.w <- strwidth(dp.labels, cex = dp.legend.cex)
	  dp.h <- strheight(dp.labels, cex = dp.legend.cex)
	  if (logic){
		 text(par("usr")[1], par("usr")[4] - 0.5*max(dp.h[1:2]) - 0.1*par("cxy")[2], parse(text = paste("list(", 
			deparse(dp.labels[[1]], width.cutoff = 500L), ", ~~", deparse(dp.labels[[2]], width.cutoff = 500L), ")", sep="")), 
			adj = 0, cex = dp.legend.cex)
	  	 text(par("usr")[1], par("usr")[4] - max(dp.h[1:2]) - 0.5*dp.h[3] - 0.2*par("cxy")[2], dp.labels[3], adj = 0, cex = dp.legend.cex)
		 y3 <- y3 - max(dp.h[1:2]) - dp.h[3] - 0.3*par("cxy")[2]		 
	  }
	  else{
		 text(par("usr")[1], par("usr")[4] - 0.5 * dp.h[1] - 0.1*par("cxy")[2], dp.labels[1], adj = 0, cex = dp.legend.cex)
		 text(par("usr")[1], par("usr")[4] - dp.h[1] - 0.5*dp.h[2] - 0.2*par("cxy")[2], parse(text = paste("paste(phantom(\" Diffusion process: \"), ~~", 
			deparse(dp.labels[[2]], width.cutoff = 500L), ")", sep="")), adj = 0, cex = dp.legend.cex)
		 text(par("usr")[1], par("usr")[4] - sum(dp.h[1:2]) - 0.5*dp.h[3] - 0.3*par("cxy")[2], dp.labels[3], adj = 0, cex = dp.legend.cex)
		 y3 <- y3 - sum(dp.h[1:2]) - dp.h[3] - 0.4*par("cxy")[2]		 
	  }
	  abline(h = y3)          
    }

    ticks <- axTicks(2)
    axis(2, at = ticks[ticks <= y3], las = 1, mgp = c(3, 0.35 + 0.4/sqrt(par("cex")), 0), tcl = -0.35)
    axis(1, mgp = c(3, 0.35, 0), tcl = -0.35)
    title(main = "Approximate First-Passage-Time Density Function Plot", line = 0.85)
    if (instants) title(xlab = "t", line = 4) else title(xlab = "t", line = 1.6)
    if (ylab) title(ylab = parse(text = "g[1](t)"), line = 3)
        
    if (any(growth.points, instants)) {
        Y <- attr(x, "summary.fptl")
        if (growth.points) {
            i <- which(Y[, 1] >= x$x[1])
            if (length(i) > 0) {
                segments(Y[i, 1], par("usr")[3], Y[i, 1], y3, col = "darkgray", lwd = 1)
                text(Y[i, 1] - 0.1 * par("cxy")[1], y3 - 0.5*strheight(expression(t), cex = sqrt(par("cex"))) - 0.1*par("cxy")[2], 
				parse(text = paste("t[~ ", i, "]", sep = "")), adj = 1, col = "darkgray")
            }
        }
        if (instants) {
            x.t <- matrix(Y[, c(2, 3, 5)], ncol = 3)
            segments(x.t, par("usr")[3], x.t, y3, lty = 8, lwd = 1)
            text(Y[, 2] - 0.1 * par("cxy")[1], - 0.65*strheight(expression(t[1]), cex = sqrt(par("cex"))) - 0.1*par("cxy")[2],
			parse(text = paste("t[~ ", 1:nrow(x.t), "]^{~ ", shQuote("*"), "}", sep = "")), adj = 1)
            mtext(parse(text = paste("t[list(~ max,", 1:nrow(Y), ")]^{~ ", shQuote("-"), "}", sep = "")), side = 1, line = 1.6, at = Y[, 3], 
			adj = 0, ...)
            mtext(parse(text = paste("t[list(~ max,", 1:nrow(Y), ")]^{~ ", shQuote("+"), "}", sep = "")), side = 1, line = 3, at = Y[, 5], 
			adj = 0, ...)
        }
    }
}
