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
      
    par(cex = 1, ps = 9)    
    par(mar = c(3 + 2 * instants, 3.25 + ylab, 2.5, 1) + 0.1, ...)    
    pin.height <- par("pin")[2] - sum(par("mai")[c(1,3)]*(par("cex")-1))

    if (dp.legend) {	  
        args <- as.list(attr(x, "Call"))
	  if (is.element("env", names(args))){
		if (is.call(args$env)) args <- c(args[3:6], unlist(as.list(args$env)[-1]))
		else if (length(args$env) > 0) args <- c(args[3:6], unlist(args$env)) else args <- args[3:6]
	  }
	  else args <- args[3:6]
        dp.labels <- vector("expression", 3)
        dp.labels[[1]] <- substitute(paste(" Diffusion process: ", 
            list(group("{", list(X(t), paste(t, "  in  ", group("[", 
                list(t0, T), "]"))), "}"), ~~P(X(t0) == x0) == 
                1)), args)
	  
	  logic <- unlist(lapply(args, function(x) if (is.call(x)) return(TRUE) else if (is.numeric(x)) any(x < 0, format(x) != format(x, scientific=FALSE)) else return(FALSE))) 		
	  if (any(logic)) args[logic] <- lapply(args[logic], function(x) as.call(parse(text=paste("(",deparse(x),")",sep="")))[[1]])				
	  dp.labels[[2]] <- substitute(list(paste(A[1](x, t) == m, " "), 
            ~~A[2](x, t) == v), list(m = eval(parse(text = paste("substitute(", 
            gsub("*", "%.%", attr(x, "dp")$mean, fixed=TRUE), ", args)", sep = ""))), v = eval(parse(text = paste("substitute(", 
            gsub("*", "%.%", attr(x, "dp")$var, fixed=TRUE), ", args)", sep = "")))))

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

    h <- pin.height - dp.height - gp.h - ti.h
    A <- - ti.h/h
    B <- 1 + (dp.height + gp.h)/h
    y2 <- ((B-A) + 1.08*(A+B))/2.16
    y1 <- A + B - y2
    plot(x$x, x$y, xlab = "", ylab = "", type = "l", las = 1, ylim = c(y1, y2), axes = FALSE)
    box()
    axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), las = 1, mgp = c(3, 0.35 + 0.4/sqrt(par("cex")), 0), tcl = -0.35)
    axis(1, mgp = c(3, 0.35, 0), tcl = -0.35)
    title(main = "First-Passage-Time Location Function Plot", line = 0.85)
    if (instants) title(xlab = "t", line = 4) else title(xlab = "t", line = 1.6)
    if (ylab) title(ylab = "FPTL(t)", line = 3)

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
    
    if (growth.points) {		
         segments(sfptl[, 1], par("usr")[3], sfptl[, 1], y3, col = "darkgray", lwd = 1)
	   text(sfptl[, 1] - 0.1 * par("cxy")[1], y3 - 0.5*strheight(expression(t), cex = sqrt(par("cex"))) - 0.1*par("cxy")[2], 
			parse(text = paste("t[~ ", 1:nrow(sfptl), "]", sep = "")), adj = 1, col = "darkgray")
    }

    if (instants) {
         x.t <- matrix(sfptl[, c(2, 3, 5)], ncol = 3)
         y <- matrix(attr(sfptl, "FPTLValues")[, c(2, 3, 5)], ncol = 3)
         points(x.t, y, type = "p", cex = min(1, 1/par("cex")))
         segments(x.t, par("usr")[3], x.t, y, lty = 8, lwd = 1)
         segments(par("usr")[1], y, x.t, y, lty = 8, lwd = 1)
	   text(sfptl[, 2] - 0.1 * par("cxy")[1], - 0.65*strheight(expression(t[1]), cex = sqrt(par("cex"))) - 0.1*par("cxy")[2], 
			parse(text = paste("t[~ ", 1:nrow(x.t), "]^{~ ", shQuote("*"), "}", sep = "")), adj = 1)
	   mtext(parse(text = paste("t[list(~ max,", 1:nrow(sfptl), ")]^{~ ", shQuote("-"), "}", sep = "")), side = 1, line = 1.6, at = sfptl[, 3], 
			adj = 0, ...)
         mtext(parse(text = paste("t[list(~ max,", 1:nrow(sfptl), ")]^{~ ", shQuote("+"), "}", sep = "")), side = 1, line = 3, at = sfptl[, 5], 
			adj = 0, ...)
    }    
}
