is.summary.fptl <-
function (obj) 
{
    if (inherits(obj, "summary.fptl") & is.matrix(obj) & is.numeric(obj) & 
        identical(names(attributes(obj)), c("dim", "Call", "FPTLValues", "FPTLCall", "dp", "class"))) 
        if ((ncol(obj) == 5) & is.call(attr(obj, "Call")) &  
            is.matrix(attr(obj, "FPTLValues")) & is.numeric(attr(obj, 
            "FPTLValues")) & is.call(attr(obj, "FPTLCall")) & 
            is.diffproc(attr(obj, "dp"))){
			args <- as.list(attr(obj, "Call"))[-1]
			if (is.element("zeroSlope", names(args))) zeroSlope <- eval(args$zeroSlope) else zeroSlope <- 0.01
			if (is.element("p0.tol", names(args))) p0.tol <- eval(args$p0.tol) else p0.tol <- 8
			if (is.element("k", names(args))) k <- eval(args$k) else k <- 3
			if (is.numeric(zeroSlope) & is.numeric(p0.tol) & is.numeric(k))
            		if ((zeroSlope <= 0.01) & (p0.tol >= 8) & (k >= 1.5) & all(diff(c(t(obj))) >= 0)) {
                			args <- as.list(attr(obj, "FPTLCall"))[-1]
                			if (all(is.element(c("dp", "t0", "T", "x0", "S"), 
                  		names(args)))) 
                  		return(all(obj >= eval(args$t0)) & all(obj <= eval(args$T)))
				}
            }
	
    return(FALSE)
}

