# 直接法

candidate <- 1/pi
M <- pi/2
n <- 100000

theta <- runif(n, min=0, max=pi)
u <- runif(n, min=0, max=1)

condition <- (u < sin(theta))

hist <- (theta[condition])

hist(hist, breaks=seq(0, pi, by=pi/20), 
main="Histgram of data", xlab="theta", ylab="density", freq=FALSE)
curve(sin(x)/2, 0, pi, add=TRUE) 
abline(h=1/2, col="red") 
abline(h=1/pi, col="blue") 