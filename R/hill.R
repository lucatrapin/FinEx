hill <- function(x, k){

  n <- length(x)
  q <- round((1-k)*n)
  large <- sort(x,decreasing=T)
  klarge <- large[q]
  xlarge <- large[1:(q-1)]
  xi <- mean(log(xlarge/klarge))

  return(list(xi=xi, zk=klarge))

}
