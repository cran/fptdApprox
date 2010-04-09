is.fpt.density <-
function (obj) 
{
    if (inherits(obj, "fpt.density") & is.list(obj) & (length(obj) == 
        2)) 
        if (all(unlist(lapply(obj, is.numeric))) & identical(names(obj), 
            c("x", "y")) & identical(names(attributes(obj)), 
            c("names", "Call", "Steps", "cumIntegral", "skips", 
                "CPUTime", "summary.fptl", "class"))) 
            if ((length(obj$x) == length(obj$y)) & is.call(attr(obj, 
                "Call")) & is.summary.fptl(attr(obj, "summary.fptl")) & 
                is.matrix(attr(obj, "Steps")) & is.numeric(attr(obj, 
                "Steps")) & is.numeric(attr(obj, "cumIntegral")) & 
                is.logical(attr(obj, "skips")) & is.matrix(attr(obj, 
                "CPUTime")) & is.numeric(attr(obj, "CPUTime"))) 
                if (all(diff(attr(obj, "cumIntegral")) >= 0) & 
                  (length(attr(obj, "cumIntegral")) == length(attr(obj, 
                    "skips"))) & (length(attr(obj, "cumIntegral")) == 
                  nrow(attr(obj, "CPUTime"))) & (nrow(attr(obj, 
                  "Steps")) >= nrow(attr(obj, "CPUTime"))) & 
                  (ncol(attr(obj, "Steps")) == 4) & (ncol(attr(obj, 
                  "CPUTime")) == 2)) {
                  args <- as.list(attr(obj, "Call"))[-1]
                  if (all(is.element(names(args), c("sfptl", 
                    "variableStep", "from.t0", "to.T", "skip", 
                    "n", "p", "tol")))) {
                    if (is.element("variableStep", names(args))) 
                      variableStep <- as.logical(as.character(args$variableStep))
                    else variableStep <- TRUE
                    if (is.element("from.t0", names(args))) 
                      from.t0 <- as.logical(as.character(args$from.t0))
                    else from.t0 <- FALSE
                    if (is.element("to.T", names(args))) 
                      to.T <- as.logical(as.character(args$to.T))
                    else to.T <- FALSE
                    if (is.element("skip", names(args))) 
                      skip <- as.logical(as.character(args$skip))
                    else skip <- TRUE
                    if (is.element("n", names(args))) 
                      n <- args$n
                    else n <- 250
                    if (is.element("p", names(args))) 
                      p <- args$p
                    else p <- 0.2
                    if (is.element("tol", names(args))) 
                      tol <- args$tol
                    else tol <- 1e-05
                    h <- attr(obj, "Steps")[, 3]
                    h[attr(obj, "skips")] <- attr(obj, "Steps")[attr(obj, 
                      "skips"), 2] - attr(obj, "Steps")[attr(obj, 
                      "skips"), 1]
                    h <- h[1:length(attr(obj, "cumIntegral"))]
                    lowerend <- attr(obj, "Steps")[1:length(attr(obj, 
                      "cumIntegral")), 1]
                    x0 <- lowerend[1]
                    lowerend[1] <- x0 + h[1]
                    return(all(obj$y >= 0) & identical(attr(obj, 
                      "Steps"), Integration.Steps(attr(obj, "summary.fptl"), 
                      variableStep, from.t0, to.T, n, p)) & identical(obj$x, 
                      c(x0, lowerend[1], unlist(mapply(seq, lowerend + 
                        h, attr(obj, "Steps")[1:length(attr(obj, 
                        "cumIntegral")), 2], by = h)))))
                  }
                }
    return(FALSE)
}

