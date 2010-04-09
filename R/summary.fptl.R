summary.fptl <-
function (object, zeroSlope = 0.01, p0.tol = 8, ...) 
{
    if (!is.fptl(object)) 
        stop(paste(sQuote("object"), "is not of class", shQuote("fptl")))
    if (zeroSlope > 0.01) {
        warning("zeroSlope is too big, default value has been used", 
            call. = FALSE)
        zeroSlope <- 0.01
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
        Y <- mapply(function(i, j, z) z[i:j], G[, 1], G[, 2], 
            MoreArgs = list(z = object$y), SIMPLIFY = FALSE)
        indexes[, 2] <- mapply(function(x, y) which(x >= y)[1], 
            Y, probs[, 2]) + G[, 1] - 1
        indexes[, 3] <- mapply(function(x, y) which(x >= y)[1], 
            Y, probs[, 3]) + G[, 1] - 1
        Y <- mapply(function(i, j, z) z[i:j], G[, 2], c(G[, 1][-1], 
            length(object$x)), MoreArgs = list(z = object$y), 
            SIMPLIFY = FALSE)
        indexes[, 5] <- mapply(function(x, y) if (any(x < y)) 
            return(which(x < y)[1] - 1)
        else return(which.min(x)), Y, probs[, 5]) + G[, 2] - 
            1
        probs[, 2] <- object$y[indexes[, 2]]
        probs[, 3] <- object$y[indexes[, 3]]
        probs[, 5] <- object$y[indexes[, 5]]
        out <- matrix((object$x)[indexes], ncol = 5)
        attr(out, "Call") <- match.call()
        attr(out, "zeroSlope") <- zeroSlope
        attr(out, "p0.tol") <- p0.tol
        attr(out, "FPTLValues") <- probs
        attr(out, "FPTLCall") <- attr(object, "Call")
        attr(out, "dp") <- attr(object, "dp")
        class(out) <- c("summary.fptl", "matrix")
        return(out)
    }
}

