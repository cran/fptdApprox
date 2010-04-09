Approx.fpt.density <-
function (sfptl, variableStep = TRUE, from.t0 = FALSE, to.T = FALSE, 
    skip = TRUE, n = 250, p = 0.2, tol = 1e-05) 
{
    if (!is.summary.fptl(sfptl)) 
        stop(paste(sQuote("sfptl"), " object is not of class ", 
            shQuote("summary.fptl"), ".", sep = ""))
    L <- Integration.Steps(sfptl, variableStep, from.t0, to.T, 
        n, p)
    g <- alist(x = , y = )
    attr(g, "Call") <- match.call()
    attr(g, "Call") <- call("Approx.density.fpt", sfptl = substitute(sfptl), 
        variableStep = variableStep, from.t0 = from.t0, to.T = to.T, 
        skip = skip, n = n, p = p, tol = tol)
    attr(g, "Steps") <- L
    args <- as.list(attr(sfptl, "FPTLCall"))
    args <- c(args[3:6], unlist(args[[7]]))
    rho <- sign(eval(parse(text = args$S), list(t = args$t0)) - 
        args$x0)
    if (rho == 0) 
        stop("S(t0) = x0")
    if (inherits(try(parse(text = args$S), silent = TRUE), "try-error")) 
        stop(paste("the mathematical expression of the boundary shows syntax errors.", 
            sep = ""))
    if (inherits(try(D(parse(text = args$S), "t"), silent = TRUE), 
        "try-error")) 
        stop("R can not compute the symbolic derivative with respect to 't' of the mathematical expression of the boundary")
    DB <- deriv(parse(text = args$S), "t", function.arg = "t")
    A1 <- function(x, t) NULL
    body(A1) <- eval(parse(text = paste("substitute(", attr(sfptl, 
        "dp")$mean, ", args)", sep = "")))
    DA2 <- deriv(eval(parse(text = paste("substitute(", attr(sfptl, 
        "dp")$var, ", args)", sep = ""))), "x", function.arg = c("x", 
        "t"))
    Df <- deriv(eval(parse(text = paste("substitute(", attr(sfptl, 
        "dp")$tpdf, ", args)", sep = ""))), "x", function.arg = c("x", 
        "t", "y", "s"))
    g0 <- L[1, 1]
    h <- L[1, 3]
    g$x <- L[1, 1] + h
    b <- DB(g$x)
    a2 <- DA2(b, g$x)
    ff <- Df(b, g$x, args$x0, args$t0)
    g$y <- max(0, -rho * (ff * (attr(b, "gradient") - A1(b, g$x) + 
        0.75 * attr(a2, "gradient")) + attr(ff, "gradient") * 
        a2))
    F <- h * g$y/2
    now <- matrix(, nrow = nrow(L) + 1, ncol = 2)
    now[1, ] <- proc.time()[1:2]
    Integral <- numeric(nrow(L))
    jumps <- logical(nrow(L))
    L[1, 1] <- g$x
    for (i in 1:nrow(L)) {
        if (skip & L[i, 4] & (g$y[length(g$y)] < 10^-10)) {
            correction <- h[length(h)] * g$y[length(g$y)]/2
            Integral[i - 1] <- Integral[i - 1] - correction
            F <- F - correction
            g$y[length(g$y)] <- 0
            h <- c(h, L[i, 2] - L[i, 1])
            g$x <- c(g$x, L[i, 2])
            b <- c(b, DB(L[i, 2]))
            g$y <- c(g$y, 0)
            now[i + 1, ] <- proc.time()[1:2]
            Integral[i] <- 0
            jumps[i] <- TRUE
        }
        else {
            u2 <- seq(L[i, 1] + L[i, 3], L[i, 2], by = L[i, 3])
            h2 <- rep(L[i, 3], length(u2))
            b2 <- DB(u2)
            a1 <- A1(b2, u2)
            a2 <- DA2(b2, u2)
            a <- as.vector(attr(b2, "gradient")) - a1 + 0.75 * 
                as.vector(attr(a2, "gradient"))
            f0 <- Df(b2, u2, args$x0, args$t0)
            if (length(a) == 1) 
                a <- rep(a, length(u2))
            if (length(a2) == 1) 
                a2 <- rep(a2, length(u2))
            g$y <- c(g$y, -rho * (f0 * a + as.vector(attr(f0, 
                "gradient")) * a2))
            n <- length(g$x)
            h <- c(h, h2)
            g$x <- c(g$x, u2)
            if (length(b2) == 1) 
                b2 <- rep(b2, length(u2))
            b <- c(b, b2)
            for (k in (1:length(u2))) {
                index <- 1:(n + k - 1)
                f1 <- Df(b2[k], u2[k], b[index], g$x[index])
                if (length(f1) == 1) {
                  Df1 <- attr(f1, "gradient")
                  f1 <- rep(f1, n + k - 1)
                  attr(f1, "gradient") <- Df1
                }
                g$y[n + k] <- max(0, g$y[n + k] + rho * sum(h[index] * 
                  g$y[index] * (f1 * a[k] + as.vector(attr(f1, 
                  "gradient")) * a2[k])))
            }
            Integral[i] <- sum(h2 * (g$y[n + (1:length(u2))] + 
                g$y[n - 1 + (1:length(u2))]))/2
            F <- F + Integral[i]
            now[i + 1, ] <- proc.time()[1:2]
            if ((!to.T) & (!L[i, 4]) & (F >= (1 - tol))) 
                break
        }
    }
    g$x <- c(g0, g$x)
    g$y <- c(0, g$y)
    if (i < nrow(L)) {
        Integral <- Integral[1:i]
        jumps <- jumps[1:i]
        now <- now[1:(i + 1), ]
    }
    CPUTime <- apply(now, 2, diff)
    if (!is.matrix(CPUTime)) 
        CPUTime <- matrix(CPUTime, ncol = 2)
    attr(g, "cumIntegral") <- cumsum(Integral)
    attr(g, "skips") <- jumps
    attr(g, "CPUTime") <- CPUTime
    attr(g, "summary.fptl") <- sfptl
    class(g) <- c("fpt.density", "list")
    return(g)
}

