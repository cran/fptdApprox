summary.fptl <-
function (object, zeroSlope = 0.01, p0.tol = 8, k = 3, ...) 
{
    if (!is.fptl(object)) 
        stop(paste(sQuote("object"), "is not of class", shQuote("fptl")))

    Call <- match.call()

    if (zeroSlope > 0.01) {
        warning("zeroSlope is too big, default value has been used", 
            call. = FALSE)	  
        zeroSlope <- 0.01
	  Call$zeroSlope <- 0.01
    }
    if (p0.tol < 8) {
        warning("p0.tol is too small, default value has been used", 
            call. = FALSE)
        p0.tol <- 8
	  Call$p0.tol <- 8
    }
    if (k < 1.5) {
        warning("k is too small, the value 1.5 has been used", 
            call. = FALSE)
        k <- 1.5
	  Call$k <- 1.5
    }

    G <- growth.intervals(object$x, object$y, zeroSlope)
	
    if (is.null(G)) 
        stop("the FPTL function is not growing")
    else {
        probs <- matrix(, nrow = nrow(G), ncol = 5)
        probs[, 1] <- (object$y)[G[, 1]]
        probs[, 4] <- (object$y)[G[, 2]]
        p <- probs[, 4] - probs[, 1]
        probs[, 2] <- probs[, 1] + p * (10^{
            -p0.tol
        })
        probs[, 3] <- probs[, 4] * (1 - 0.05 * p)
        probs[, 5] <- 1 - (1 - probs[, 4]^2)^(0.5 * (1 + p/probs[, 
            4]))
	
        indexes <- matrix(, nrow = nrow(G), ncol = 5)
        indexes[, c(1, 4)] <- G

	  I <- mapply(seq, from = G[, 1], to = G[, 2], SIMPLIFY = FALSE)			 
	  indexes[, 2] <- mapply(function(index,p,Y) which(Y[index] >= p)[1], I, probs[, 2], MoreArgs = list(Y = object$y)) + G[, 1] - 1	  
	  indexes[, 3] <- mapply(function(index,p,Y) which(Y[index] >= p)[1], I, probs[, 3], MoreArgs = list(Y = object$y)) + G[, 1] - 1

        probs[, 2] <- object$y[indexes[, 2]]
        probs[, 3] <- object$y[indexes[, 3]]  
	  
	  I <- mapply(seq, from = G[, 2], to = c(G[, 1][-1], length(object$x)), SIMPLIFY = FALSE)	
	  x1 <- object$x[c(G[, 1][-1], length(object$x))]
	  x2 <- (object$x)[indexes[,4]] + k*((object$x)[indexes[,4]] - (object$x)[indexes[,2]])*(1 - probs[, 4])			
	  logic <- (x2 < x1)	
	  if (any(logic)) I[which(logic)] <- mapply(function(index,s,x) index[1] + which(x[index] <= s) - 1, I[which(logic)], x2[logic], MoreArgs = list(x = object$x), SIMPLIFY = FALSE)	  
	  indexes[, 5] <- mapply(function(index, p, Y) max(which(Y[index] >= p)), I, probs[, 5], MoreArgs = list(Y = object$y)) + G[, 2] - 1
	  
	  probs[, 5] <- object$y[indexes[, 5]]
	  	  
	  out <- matrix((object$x)[indexes], ncol = 5)	
        attr(out, "Call") <- Call	  
        attr(out, "FPTLValues") <- probs
        attr(out, "FPTLCall") <- attr(object, "Call")
        attr(out, "dp") <- attr(object, "dp")

	  aux <- function(x, envir) eval(parse(text=paste("substitute(", deparse(x), ", envir)", sep="")))
	  if (length(attr(out, "Call")) > 2){   	
    	  	vars <- lapply(attr(out, "Call")[-(1:2)], all.vars)    
    	  	logic <- (unlist(lapply(vars, length)) > 0)
    	  	if (any(logic)){		
			envars <- lapply(vars[logic], mget, envir=parent.frame())				
			attr(out, "Call")[-(1:2)][logic] <- mapply(aux, attr(out, "Call")[-(1:2)][logic], envars, SIMPLIFY = FALSE)		
    	  	}
    	
    		logic <- unlist(lapply(lapply(lapply(attr(out, "Call")[-(1:2)], all.names), is.element, el=c("[","[[")), any))
    		if (any(logic)) attr(out, "Call")[-(1:2)][logic] <- lapply(attr(out, "Call")[-(1:2)][logic], eval)	
	  }

        class(out) <- c("summary.fptl", "matrix")
        return(out)
    }
}

