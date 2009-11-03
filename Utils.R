ReadParameters <- function(file)
  {
    parameters<- readLines(file)
    n.par <- length(parameters)
    par.name <- NULL
    par.value <- NULL

    
    for( i in 1:n.par)
      {
        par.i <- unlist(strsplit(parameters[i], "#"))
        par.i <- unlist(strsplit(par.i[1], "  "))
        par.name=c(par.name,par.i[1])
        par.values <- list(type.convert(unlist(strsplit(par.i[2], " "))))
        par.value=c(par.value,par.values)
      }
    para <- as.list(par.value)
    names(para) <- par.name
    return(para)
  }


Add2File <- function(x, file.out, append=TRUE)
  {
    nb.col=length(x)
    write(x, file=file.out, ncolumns=nb.col, sep=",", append=append)
  }
