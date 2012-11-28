FPTL <-
function (dp, t0, T, x0, S, env = NULL, n = 10000) 
{	
    if (!is.diffproc(dp)) 
        stop(paste(sQuote("dp"), "argument is not of class ", shQuote("diffproc"), 
            ".", sep = ""))
    if (!is.numeric(t0)) stop(paste(sQuote("t0"), "argument is not of numeric type."))
    if (length(t0) > 1) stop(paste(sQuote("t0"), "argument is not of length 1."))
    if (!is.numeric(T)) stop(paste(sQuote("T"), "argument is not of numeric type."))
    if (length(T) > 1) stop(paste(sQuote("T"), "argument is not of length 1."))    
    if (!is.numeric(x0)) stop(paste(sQuote("x0"), "argument is not of numeric type."))
    if (length(x0) > 1) stop(paste(sQuote("x0"), "argument is not of length 1."))
    if (!is.character(S)) stop(paste(sQuote("S"), "argument is not of character type."))
    if (length(S) > 1) stop(paste(sQuote("S"), "argument is not of length 1."))
    if ((!is.list(env)) & (!is.null(env))) 
        stop(paste(sQuote("env"), "argument is not NULL or a list object."))
    if (length(env) == 0)
        warning("env argument is not specified. R looks for the objects not found into the temporary environment for evaluating dp object in the parent.frame() environment.")
    if (!is.numeric(n)) stop(paste(sQuote("n"), "argument is not of numeric type."))
    if (round(n) != n) stop(paste(sQuote("n"), "argument is not an integer number."))
    if (length(n) > 1) stop(paste(sQuote("n"), "argument is not of length 1."))
    if (t0 >= T) 
        stop("the final time instant is not greater than the initial time instant.")
    if (inherits(try(parse(text=S), silent = TRUE), "try-error")) 
        stop(paste("the mathematical expression of the boundary shows syntax errors.", 
            sep = ""))   
    if (inherits(try(D(parse(text=S), "t"), silent = FALSE), "try-error")) 
        stop("R can not compute the symbolic derivative with respect to 't' of the mathematical expression of the boundary")

    s0 <- eval(parse(text = S), c(env, list(t = t0, t0 = t0, x0 = x0)))	
    if (x0 == s0) 
        stop("the value of the boundary at the initial time instant is equal to the initial value of the process.")

    grid.t <- seq(t0, T, length = n)[-1]
    S.t <- eval(parse(text = S), c(env, list(t = grid.t, t0 = t0, x0 = x0)))
	
    y.t <- eval(parse(text = dp[[4]]), c(env, list(x = S.t, t = grid.t, 
        y = x0, s = t0)))
    if (x0 < s0) 
        y.t <- 1 - y.t
	
    G <- growth.intervals(grid.t, y.t)
    if (is.null(G)) 
        warning("the FPTL function is not growing.")
    else {
        index <- c(t(G))
        endlimits <- grid.t[index]
        m <- diff(endlimits)/(T - t0)
        if (G[1, 1] > 1) {
            endlimits <- c(t0, endlimits)
            m <- c((grid.t[G[1, 1]] - t0)/(5 * (T - t0)), m)
        }
        if (G[length(G)] < (n - 1)) {
            i <- which.min(y.t[(1 + G[length(G)]):(n - 1)])
            endlimits <- c(endlimits, grid.t[i + G[length(G)]])
            m <- c(m, (grid.t[i + G[length(G)]] - grid.t[G[length(G)]])/(T - 
                t0))
            j <- G[length(G)] + i
            if (j < n) {
                endlimits <- c(endlimits, T)
                m <- c(m, (T - grid.t[i + G[length(G)]])/(5 * 
                  (T - t0)))
            }
        }
        m <- trunc(m * (n - 1)/sum(m))
        m <- (m + 1 + abs(m - 1))/2
        d <- (n - 1) - sum(m)
        if (d > 0) {
            j <- rep(order(m), length.out = d)
            m[j] <- m[j] + 1
        }
        if (d < 0) {
            i <- order(m, decreasing = T)
            j <- rep(i[m[i] > d], length.out = d)
            m[j] <- m[j] - 1
        }
        d <- (n - 1) - sum(m)
        if (d != 0) 
            stop("n is too small.")
        if (length(endlimits) == 2) 
            grid.t <- seq(endlimits[1], endlimits[2], length.out = m + 
                1)[-1]
        else {
            v <- mapply(seq, endlimits[-length(endlimits)], endlimits[-1], 
                length.out = m + 1, SIMPLIFY = FALSE)
            if (is.list(v)) 
                grid.t <- unlist(lapply(v, function(l) l[-1]))
            else grid.t <- as.vector(apply(v, 2, function(l) l[-1]))
        }
        S.t <- eval(parse(text = S), c(env, list(t = grid.t, t0 = t0, x0 = x0)))
        y.t <- eval(parse(text = dp[[4]]), c(env, list(x = S.t, 
            t = grid.t, y = x0, s = t0)))
        if (x0 < s0) 
            y.t <- 1 - y.t
    }		

    fptl <- list(x = c(t0, grid.t), y = c(0, y.t))

    attr(fptl, "dp") <- dp        
    attr(fptl, "Call") <- match.call()
		    	
    aux <- function(x, envir) eval(parse(text=paste("substitute(", deparse(x), ", envir)", sep="")))
    vars <- lapply(attr(fptl, "Call")[-(1:2)], all.vars)
    logic <- (unlist(lapply(vars, length)) > 0)
    if (any(logic)){		
		envars <- lapply(vars[logic], mget, envir=parent.frame())				
		attr(fptl, "Call")[-(1:2)][logic] <- mapply(aux, attr(fptl, "Call")[-(1:2)][logic], envars, SIMPLIFY = FALSE)		
    }
    
    if (missing(n)) index <- 3:6 else index <- c(3:6, 8 - missing(env))		
    logic <- unlist(lapply(lapply(lapply(attr(fptl, "Call")[index], all.names), is.element, el=c("[","[[")), any))
    if (any(logic)) attr(fptl, "Call")[index][logic] <- lapply(attr(fptl, "Call")[index][logic], eval)	

    attr(fptl, "Call")$S <- as.character(parse(text=attr(fptl, "Call")$S))

    if (length(env) > 0){
	if (is.call(attr(fptl, "Call")$env)){
		logic <- unlist(lapply(lapply(lapply(attr(fptl, "Call")$env[-1], all.names), is.element, el=c("[","[[")), any))
		if (any(logic)) attr(fptl, "Call")$env[-1][logic] <- lapply(attr(fptl, "Call")$env[-1][logic], eval)
	}	
    }
    
    class(fptl) <- c("fptl", "list")
    return(fptl)
}
