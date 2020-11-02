# Gibbs sampling

# 事前分布のハイパーパラメータ
mu0 <- 50; sig0 <-10
r0 <- 4; s0 <- 200

iterations <- 3500
y <- c(60,100,76,35,14,25,56,30,22,50,33,76,58,67,28,48,63,53,54,33,34,38,33,57,63,84,46,28,79,89)
ybar <- mean(y)
n <- length(y)
mat <- matrix(nrow=iterations, ncol=2)

# 初期値
mu <- rnorm(1, mean=mu0, sd=sig0)
tau <- rgamma(1, r0/2, s0/2)
sig <- sqrt(1/tau) # stdであることに注意

for(t in 1:iterations){
	w <- sig0^2 / ((sig^2)/n+sig0^2)
	m <- w*ybar + (1-w)*mu0
	s <- sqrt(w*(sig^2)/n)
	mu <- rnorm(1, m, s)

	a <- r0 + n
	b <- s0 + sum((y-mu)^2)
	tau <- rgamma(1, a/2, b/2)
	sig <- sqrt(1/tau)

	mat[t,] <- c(mu, sig)
}
burn.out <- mat[501:iterations,]
op <- par(mfrow=c(2,2))
p <- 1:iterations

plot(p, mat[,1], type="l", main="mu", xlab=" iteration", ylab="mu", lty=1, lwd=1);
plot(p, mat[,2], type="l", main="sigma", xlab=" iteration", ylab="sigma", lty=1, lwd=1);
hist(burn.out[,1])
hist(burn.out[,2])
par(op)

# 数値出力
cat("mu Bayesian estimation   :")
print(mean(burn.out[,1]))
cat("mu MLE                   :")
print(ybar)
cat("sigma Bayesian estimation:")
print(mean(burn.out[,2]))
cat("variance                 :")
print(var(y))
cat("std                      :")
print(sd(y))