FPTL <-
function (dp, t0, T, x0, S, env = NULL, n = 10000) 
{
    if (!is.diffproc(dp)) 
        stop(paste(sQuote("dp"), "object is not of class ", shQuote("diffproc"), 
            ".", sep = ""))
    if ((!is.list(env)) & (!is.null(env))) 
        stop(paste(sQuote("env"), "argument is not NULL or a list object."))
    if (is.null(env)) 
        warning("env argument is not specified. R looks for the objects not found into the temporary environment \n\t\t\t\tfor evaluating dp object in the parent.frame() environment.")
    if (t0 >= T) 
        stop("the final time instant is not greater than the initial time instant.")
    if (inherits(try(parse(text = S), silent = TRUE), "try-error")) 
        stop(paste("the mathematical expression of the boundary shows syntax errors.", 
            sep = ""))
    s0 <- eval(parse(text = S), list(t = t0))
    if (x0 == s0) 
        stop("the value of the boundary at the initial time instant is equal to the initial value of the process.")
    grid.t <- seq(t0, T, length = n)[-1]
    S.t <- eval(parse(text = S), list(t = grid.t))
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
                length.out = m + 1, simplify = FALSE)
            if (is.list(v)) 
                grid.t <- unlist(lapply(v, function(l) l[-1]))
            else grid.t <- as.vector(apply(v, 2, function(l) l[-1]))
        }
        S.t <- eval(parse(text = S), list(t = grid.t))
        y.t <- eval(parse(text = dp[[4]]), c(env, list(x = S.t, 
            t = grid.t, y = x0, s = t0)))
        if (x0 < s0) 
            y.t <- 1 - y.t
    }
    fptl <- list(x = c(t0, grid.t), y = c(0, y.t))
    attr(fptl, "dp") <- dp
    attr(fptl, "Call") <- call("FPTL", dp = substitute(dp), t0 = t0, 
        T = T, x0 = x0, S = S, env = env, n = n)
    class(fptl) <- c("fptl", "list")
    return(fptl)
}

