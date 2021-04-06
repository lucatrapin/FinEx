arq.sv <- function(r, p, x){

  a <- utils::combn(seq(0.05,0.95,by=0.05),2)
  a <- t(a[,apply(a,2,sum)<1])
  a <- rbind(a,cbind(a[,2],a[,1]))
  n <- nrow(a)
  o <- stats::quantile(r, p)*(1 - apply(a,1,sum))

  tmp <- cbind(o, a, rep(0,n))
  for (i in 1:n){
    tmp[i,4] <- arq.llk(tmp[1:3], r, p, x)
  }

  sv <- tmp[which.min(tmp[,4]),1:3]

  return(sv)

}
