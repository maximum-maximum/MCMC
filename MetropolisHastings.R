# Metropolis-Hastings

# 事前分布のハイパーパラメータ
mu0 <- 50; sig0 <-10
r0 <- 4; s0 <- 200
fixedsig <- 5 

iterations <- 3500
y <- c(60,100,76,35,14,25,56,30,22,50,33,76,58,67,28,48,63,53,54,33,34,38,33,57,63,84,46,28,79,89)
ybar <- mean(y)
s <- sd(y)
n <- length(y)
mat <- matrix(nrow=iterations, ncol=2)

# 初期値
tau <- rgamma(1, r0/2, s0/2)
sig <- sqrt(1/tau) # stdであることに注意
mu <- rnorm(1, mean=mu0, sd=sig)

for(t in 1:iterations){
	mu1 <- (n*ybar + mu0) / (n+1) # fix
	sig1 <- sqrt(sig^2 / (n+1)) # 初期値sig依存
	r1 <- r0 + n # fix
	s1 <- s0 + n*s^2 + (n*(ybar-mu0)^2) / (n+1) # 初期値sig依存
	
	# about mu
	z <- rnorm(1, 0, fixedsig)
	u <- runif(1, 0, 1)
	mudash <- mu + z

	pi <-     (sig1^2)^(-1/2)*exp((-(mu-mu1)^2)/(2*sig1^2))*
			(sig^2)^(-r1/2-1)*exp(-s1/(2*sig^2))
	pidash <- (sig1^2)^(-1/2)*exp((-(mudash-mu1)^2)/(2*sig1^2))*
			(sig^2)^(-r1/2-1)*exp(-s1/(2*sig^2))
	alpha <- min(1, pidash/pi)

	if(u <= alpha){
		mu <- mudash
	}

	# about sigma
	z <- rnorm(1, 0, fixedsig)
	u <- runif(1, 0, 1)
	sigdash <- sig + z

	pi <-     (sig1^2)^(-1/2)*exp((-(mu-mu1)^2)/(2*sig1^2))*
			(sig^2)^(-r1/2-1)*exp(-s1/(2*sig^2))

	sig1 <- sqrt(sigdash^2/(n+1))
	s1 <- s0+n*s^2+(n*(ybar-mu0)^2)/(n+1)

	pidash <- (sig1^2)^(-1/2)*exp((-(mu-mu1)^2)/(2*sig1^2))*
			(sigdash^2)^(-r1/2-1)*exp(-s1/(2*sigdash^2))
	alpha <- min(1, pidash/pi)

	if(u <= alpha){
	  # sig1，s1の更新はforの最初にsig依存で行うため不要
		sig <- sigdash 
	}
	mat[t,] <- c(mu, sig)
}
init <- 501
burn.out <- mat[init:iterations,]
op <- par(mfrow=c(3,2))
p <- 1:iterations
q <- init:iterations
plot(p, mat[,1], type="l", main="mu", xlab=" iteration", ylab="mu", lty=1, lwd=1);
plot(p, mat[,2], type="l", main="sigma", xlab=" iteration", ylab="sigma", lty=1, lwd=1);
plot(q, burn.out[,1], type="l", main="mu", xlab=" iteration", ylab="mu", lty=1, lwd=1);
plot(q, burn.out[,2], type="l", main="sigma", xlab=" iteration", ylab="sigma", lty=1, lwd=1);
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