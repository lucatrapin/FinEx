link.fun.inv <- function(type){

  if (type=="identity"){
    out <- function(x){return(identity(x))}
  }

  if (type=="exp"){
    out <- function(x){return(log(x))}
  }

  if (type=="logit"){
    out <- function(x){return(log(x/(1-x)))}
  }

  return(out)

}
