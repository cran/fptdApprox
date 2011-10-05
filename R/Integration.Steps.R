Integration.Steps <-
function (sfptl, variableStep = TRUE, from.t0 = FALSE, to.T = FALSE, 
    n = 250, p = 0.2, alpha = 1) 
{
    if (!is.summary.fptl(sfptl)) 
        stop(paste(sQuote("sfptl"), " object is not of class ", 
            shQuote("summary.fptl"), ".", sep = ""))    

    if (n<=0) stop("n <= 0.")
    if (!variableStep){ 	  
        if (missing(p)){ 
		if (!missing(alpha)) cat("alpha argument is ignored.\n")
	  }
	  else{
		if (missing(alpha)) cat("p argument is ignored.\n")
		else cat("p and alpha arguments are ignored.\n")
	  }
    }

    t0 <- eval(attr(sfptl, "FPTLCall")$t0)
    T <- eval(attr(sfptl, "FPTLCall")$T)    

    h <- (sfptl[, 3] - sfptl[, 2])/n

    m <- 2 * nrow(sfptl) - 1        
    L <- matrix(, nrow = m, ncol = 4)
    L[, 4] <- FALSE
    odd <- seq(1, m, by = 2)	

    if (variableStep){	  
	  if ((p <= 0) | (p > 1)) stop("p is not in (0,1].")
	  if ((alpha <= 0) | (alpha > 1)) stop("alpha is not in (0,1].")
	  
	  cotes <- c(n = 50, p=0.1, alpha = 0.75)
	  logic <- (c(n,p,alpha) < cotes)
	  if (any(logic)){
		text <- paste(names(cotes)[logic], cotes[logic], sep=" >= ")		
		if (length(text) > 1) text <- paste(paste(text[-length(text)], collapse=", "), text[length(text)], sep=" and ")		
		text <- paste(text, "is recommended. If not, some integration steps can be too large.")
		cat(text, "\n")	
		repeat{ 
			answer <- readline("Do you want to continue? (y/n) ")
			if ((answer == "n") | (answer == "y")) break
		}   	
		if (answer == "n") stop("Approximation process stopped by the usuary.\n\n") else cat("\n")       	  		
	  }	  
    	  	       	               
        L[odd, 1] <- sfptl[, 2]
        L[odd, 2] <- sfptl[, 2] + h * trunc((sfptl[, 5] - sfptl[, 2])/h)
        L[odd, 3] <- h

        if (m > 1) {
		even <- seq(2, m, by = 2)
            L[even, 1] <- L[odd, 2][-nrow(sfptl)]
            L[even, 2] <- L[odd, 1][-1]
		s <- (L[even, 2] - L[even, 1])/L[even-1, 3]					                       	
		L[even, 3] <- (L[even, 2] - L[even, 1])/pmax.int(1, trunc(ifelse(s <= n, s*p, n*p*(s/n)^alpha)))		
            L[even, 4] <- TRUE		
        }
        
        if (from.t0 & (t0 < L[1, 1])){
			s <- (L[1, 1] - t0)/L[1, 3]													         
			L <- rbind(c(t0, L[1, 1], (L[1, 1] - t0)/max(1, trunc(ifelse(s <= n, s*p, n*p*(s/n)^alpha))), FALSE), L)
			
	  }
        if (to.T & (sfptl[nrow(sfptl), 5] < T)){
			s <- (T - L[nrow(L), 2])/L[nrow(L), 3]						      				
			L <- rbind(L, c(L[nrow(L), 2], T, (T - L[nrow(L), 2])/max(1, trunc(ifelse(s < n, p*s, n*p*(s/n)^alpha))), FALSE))			
	  }            	
    }
    else {
        h <- min(h)
        if (from.t0) 
            origin <- t0
        else origin <- sfptl[1, 2]
        L[odd, 1] <- origin + h * trunc((sfptl[, 2] - origin)/h)
        L[odd, 2] <- origin + h * trunc((sfptl[, 5] - origin)/h)

        if (m > 1) {
		even <- seq(2, m, by = 2)
            L[even, 1] <- L[odd, 2][-nrow(sfptl)]
            L[even, 2] <- L[odd, 1][-1]
            L[even, 4] <- TRUE
        }

        L[, 3] <- h

        if (from.t0 & (t0 < L[1, 1])) 
            L <- rbind(c(t0, L[1, 1], h, FALSE), L)
        if (to.T & (sfptl[nrow(sfptl), 5] < T)) 
            L <- rbind(L, c(L[nrow(L), 2], L[nrow(L), 2] + h * 
                trunc((T - L[nrow(L), 2])/h), h, FALSE))
    }
    index <- which(L[, 1] == L[, 2])
    if (length(index)>0) L <- L[- index, ]    	
    return(L)
}

