link.fun <- function(type){

  if (type=="identity"){
    out <- function(x){return(identity(x))}
  }

  if (type=="exp"){
    out <- function(x){return(exp(x))}
  }

  if (type=="logit"){
    out <- function(x){return(exp(x)/(1+exp(x)))}
  }

  return(out)

}
