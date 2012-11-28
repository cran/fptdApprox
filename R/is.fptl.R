is.fptl <-
function (obj) 
{
    if (inherits(obj, "fptl") & is.list(obj) & (length(obj) == 
        2)) 
        if (all(unlist(lapply(obj, is.numeric))) & identical(names(obj), 
            c("x", "y")) & identical(names(attributes(obj)), 
            c("names", "dp", "Call", "class"))) 
            if ((length(obj$x) == length(obj$y)) & is.diffproc(attr(obj, 
                "dp")) & is.call(attr(obj, "Call"))) {			
                args <- as.list(attr(obj, "Call"))
                if (any(!is.element(c("dp", "t0", "T", "x0", "S"), names(args)))) return(FALSE)
		    args <- lapply(args[-(1:2)], eval)		
		    if (any(!(unlist(lapply(args[1:3], is.numeric))))) return(FALSE)
		    if (!is.character(args$S)) return(FALSE)
		    if (any(unlist(lapply(args[1:4], length)) > 1)) return(FALSE)
		    if (inherits(try(parse(text=args$S), silent = TRUE), "try-error")) return(FALSE)
  		    if (inherits(try(D(args$S, "t"), silent = TRUE), "try-error")) return(FALSE)		
		    if (is.element("env", names(args))){
			env <- eval(attr(obj, "Call")$env)
			if (is.null(env) | is.list(env)){
				if (length(env) > 0){								
					if (all(unlist(lapply(env, is.numeric)))){									
						if (any(unlist(lapply(env, length)) > 1)) return(FALSE)
					}
					else return(FALSE)
				}
			}
			else return(FALSE)
		    }
		    if (is.element("n", names(args))) 
                  n <- eval(args$n)
                else n <- 10000
                return((length(obj$x) == n) & all(obj$x >= 
                  	args$t0) & all(obj$x <= args$T) & all(obj$y >= 
                    	0) & all(obj$y <= 1))		
            }
    return(FALSE)
}
