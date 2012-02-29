TUWmodel <- function (prec, airt, ep, area, param=c(1.2,1.2,2,-2,0,0.9,100,3.3,0.5,9,105,50,2,10,26.5), incon=c(50,0,2.5,2.5), itsteps=NULL) {
 nzones <- ifelse(is.vector(prec), 1, dim(prec)[2])
 itsteps <- ifelse(is.null(itsteps), length(prec)/nzones, itsteps)
 if (nzones == 1) {
  input <- as.matrix(t(cbind(prec, airt, ep, ep*0)))
 } else {
  input <- matrix(NA, ncol=dim(prec)[1], nrow=4*nzones)
  for (i in 1:nzones) {
   input[4*(i-1)+c(1:4),] <- t(cbind(prec[,i], airt[,i], ep[,i], ep[,i]*0))
  }
 }
 storage.mode(input) <- "double"
 output <- array(777, dim=c(nzones, 20, itsteps))
 storage.mode(output) <- "double"
 dummy <- .Fortran("hbvmodel", itsteps=as.integer(itsteps), nzones=as.integer(nzones), area=as.single(area), param=as.single(param), incon=as.single(incon),
                    input, output=output, PACKAGE="TUWmodel")
 names(dummy$param) <- c("SCF","DDF","Tr","Ts","Tm","LPrat","FC","BETA","k0","k1","k2","lsuz","cperc","bmax","croute")
 names(dummy$incon) <- c("SSM0","SWE0","SUZ0","SLZ0")
 dummy$qzones <- t(dummy$output[,1,])
 if (nzones > 1) {dummy$q <- apply(dummy$qzones,1,sum)} else {dummy$q <- dummy$qzones}
 dummy$swe <- t(dummy$output[,2,])
 dummy$q0 <- t(dummy$output[,7,])
 dummy$q1 <- t(dummy$output[,8,])
 dummy$q2 <- t(dummy$output[,9,])
 dummy$rmoist <- t(dummy$output[,3,])
 dummy$rain <- t(dummy$output[,4,])
 dummy$snow <- t(dummy$output[,5,])
 dummy$eta <- t(dummy$output[,10,])

 return(dummy)
}







