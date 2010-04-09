Integration.Steps <-
function (sfptl, variableStep = TRUE, from.t0 = FALSE, to.T = FALSE, 
    n = 250, p = 0.2) 
{
    if (!is.summary.fptl(sfptl)) 
        stop(paste(sQuote("sfptl"), " object is not of class ", 
            shQuote("summary.fptl"), ".", sep = ""))
    args <- as.list(attr(sfptl, "FPTLCall"))
    h <- (sfptl[, 3] - sfptl[, 2])/n
    m <- 2 * nrow(sfptl) - 1
    odd <- seq(1, m, by = 2)
    if (m > 1) 
        even <- seq(2, m, by = 2)
    L <- matrix(, nrow = m, ncol = 4)
    L[, 4] <- FALSE
    if (variableStep) {
        origin <- sfptl[, 2]
        L[odd, 1] <- sfptl[, 2]
        L[odd, 2] <- origin + h * trunc((sfptl[, 5] - origin)/h)
        L[odd, 3] <- h
        if (m > 1) {
            L[even, 1] <- L[odd, 2][-nrow(sfptl)]
            L[even, 2] <- L[odd, 1][-1]
            m <- round(n * p)
            L[even, 3] <- (L[even, 2] - L[even, 1])/m
            L[even, 4] <- TRUE
        }
        else m <- round(n * p)
        if (from.t0 & (args$t0 < L[1, 1])) 
            L <- rbind(c(args$t0, L[1, 1], (L[1, 1] - args$t0)/m, 
                FALSE), L)
        if (to.T & (sfptl[nrow(sfptl), 5] < args$T)) 
            L <- rbind(L, c(L[nrow(L), 2], args$T, (args$T - 
                L[nrow(L), 2])/m, FALSE))
    }
    else {
        h <- min(h)
        if (from.t0) 
            origin <- args$t0
        else origin <- sfptl[1, 2]
        L[odd, 1] <- origin + h * trunc((sfptl[, 2] - origin)/h)
        L[odd, 2] <- origin + h * trunc((sfptl[, 5] - origin)/h)
        if (m > 1) {
            L[even, 1] <- L[odd, 2][-nrow(sfptl)]
            L[even, 2] <- L[odd, 1][-1]
            L[even, 4] <- TRUE
        }
        L[, 3] <- h
        if (from.t0 & (args$t0 < L[1, 1])) 
            L <- rbind(c(args$t0, L[1, 1], h, FALSE), L)
        if (to.T & (sfptl[nrow(sfptl), 5] < args$T)) 
            L <- rbind(L, c(L[nrow(L), 2], L[nrow(L), 2] + h * 
                trunc((args$T - L[nrow(L), 2])/h), h, FALSE))
    }
    return(L)
}

