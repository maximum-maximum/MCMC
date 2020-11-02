# 間接法

n <- 100000
theta <- runif(n, min=0, max=2)
weight <- theta/sum(theta)
theta.star <- sample(theta, n, replace=TRUE, prob=weight)

hist(theta.star, breaks=seq(0, 2, 1/10),
main="Histgram of data", xlab="theta", ylab="density", freq=FALSE)
abline(h=1/2, col="blue") 
abline(0, 1/2, col="red") 
