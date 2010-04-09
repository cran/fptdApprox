is.summary.fptl <-
function (obj) 
{
    if (inherits(obj, "summary.fptl") & is.matrix(obj) & is.numeric(obj) & 
        identical(names(attributes(obj)), c("dim", "Call", "zeroSlope", 
            "p0.tol", "FPTLValues", "FPTLCall", "dp", "class"))) 
        if ((ncol(obj) == 5) & is.call(attr(obj, "Call")) & is.numeric(attr(obj, 
            "zeroSlope")) & is.numeric(attr(obj, "p0.tol")) & 
            is.matrix(attr(obj, "FPTLValues")) & is.numeric(attr(obj, 
            "FPTLValues")) & is.call(attr(obj, "FPTLCall")) & 
            is.diffproc(attr(obj, "dp"))) 
            if ((attr(obj, "zeroSlope") <= 0.01) & (attr(obj, 
                "p0.tol") > 0) & all(diff(c(t(obj))) >= 0) & 
                all(apply(matrix(attr(obj, "FPTLValues")[, 1:4], 
                  ncol = 4), 1, diff) >= 0) & all(attr(obj, "FPTLValues")[, 
                4] >= attr(obj, "FPTLValues")[, 5])) {
                args <- as.list(attr(obj, "FPTLCall"))[-1]
                if (all(is.element(c("dp", "t0", "T", "x0", "S"), 
                  names(args)))) 
                  return(all(obj >= args$t0) & all(obj <= args$T))
            }
    return(FALSE)
}

