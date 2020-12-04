### R code to estimate R-squared from the c-statistic
### useful for pre-specifying R2 for calculating sample size
### Riley RD, van Calster B, Collins GS. A note on estimating the Cox-Snell R2
### from a reported C-statistics (AUROC) to inform sample size calculations
### for developing a prediction model with a binary outcome
###
### written by Gary Collins (07-November-2020)

approximate_R2 <- function(auc, prev, n = 1000000){
  
  # define mu as a function of the C-statistic
  mu   <- sqrt(2) * qnorm(auc)
  
  # simulate large sample linear prediction based on two normals 
  # for non-eventsN(0, 1), events and N(mu, 1)
  
  LP <- c(rnorm(prev*n,  mean=0, sd=1), rnorm((1-prev)*n, mean=mu, sd=1))
  y <- c(rep(0, prev*n), rep(1, (1-prev)*n))
  
  # Fit a logistic regression with LP as covariate;
  # this is essentially a calibration model, and the intercept and    
  # slope estimate will ensure the outcome proportion is accounted 
  # for, without changing C-statistic
  
  fit <- lrm(y~LP)
  
  max_R2 <- function(prev){
    1-(prev^prev*(1-prev)^(1-prev))^2
  }
  return(list(R2.nagelkerke = as.numeric(fit$stats['R2']), 
              R2.coxsnell = as.numeric(fit$stats['R2']) * max_R2(prev)))   
}

set.seed(1234)
approximate_R2(auc = 0.81, prev = 0.77, n=1000000)

## let's see how well the approximation works
set.seed(12491)

generate_data <- function(NN, mean.lp = 0, sd.lp = 1){
  lp <- rnorm(NN, mean.lp, sd.lp)
  y <- runif(NN) < 1 / (1 + exp(lp))
  
  DATA   <- data.frame(lp)
  DATA$y <- y * 1
  DATA
}

N.SIM <- 100

#c.DEV      <- vector(mode = 'numeric', length = N.SIM)
#R2.DEV     <- vector(mode = 'numeric', length = N.SIM)
#prev.DEV   <- vector(mode = 'numeric', length = N.SIM)
#est.R2.DEV <- vector(mode = 'numeric', length = N.SIM)

N <- c(100, 250, 500, 1000)
c.DEV      <- matrix(ncol = length(N), nrow = N.SIM)
R2.DEV     <- matrix(ncol = length(N), nrow = N.SIM)
prev.DEV   <- matrix(ncol = length(N), nrow = N.SIM)
est.R2.DEV <- matrix(ncol = length(N), nrow = N.SIM)

max_R2 <- function(prev){
  1-(prev^prev*(1-prev)^(1-prev))^2
}

set.seed(2345)
for(j in 1: length(N)){
  for(i in 1:N.SIM){
    mean.lp       <- runif(1, 0, 5)
    sd.lp         <- runif(1, 0.5, 3)
    DATA          <- generate_data(N[j], mean.lp = mean.lp, sd.lp = sd.lp)
    if(length(unique(DATA$y))>1){
      prev.DEV[i, j]   <- as.numeric(prop.table(table(DATA$y))[2])  # outcome fraction
  
    # fit model
      fit.full         <- lrm(y~., data = DATA, x = T, y = T)
      c.DEV[i, j]      <- as.numeric(fit.full$stats['C']) 
      R2.DEV[i, j]     <- fit.full$stats['R2'] * max_R2(prev.DEV[i,j])
      est.R2.DEV[i, j] <- approximate_R2(auc = c.DEV[i,j], prev = prev.DEV[i,j])$R2.coxsnell # Pull out R2
    } else {
      prev.DEV[i,j]   <- NA
      c.DEV[i,j]      <- NA
      R2.DEV[i,j]     <- NA
      est.R2.DEV[i,j] <- NA
    }
  }
}

op <- par(mar = c(5,4.5,4,2) + 0.1)
par(mfrow = c(2, 2))
plot(R2.DEV[,1], est.R2.DEV[,1], 
     pch  = 20, 
     xlab = expression(paste("model ", italic(R)) ["cs"]^"2"), 
     ylab = expression(paste("approximated ", italic(R)) ["cs"]^"2"), 
     ylim = c(floor(100 * min(est.R2.DEV, R2.DEV, na.rm=T)) / 100, ceiling(100 * max(est.R2.DEV, R2.DEV, na.rm=T)) / 100), 
     xlim = c(floor(100 * min(est.R2.DEV, R2.DEV, na.rm=T)) / 100, ceiling(100 * max(est.R2.DEV, R2.DEV, na.rm=T)) / 100),
     main = "N=100")
grid()
abline(a = 0, b = 1, col = 'red')

plot(R2.DEV[,2], est.R2.DEV[,2], 
     pch  = 20, 
     xlab = expression(paste("model ", italic(R)) ["cs"]^"2"), 
     ylab = expression(paste("approximated ", italic(R)) ["cs"]^"2"), 
     ylim = c(floor(100 * min(est.R2.DEV, R2.DEV, na.rm=T)) / 100, ceiling(100 * max(est.R2.DEV, R2.DEV, na.rm=T)) / 100), 
     xlim = c(floor(100 * min(est.R2.DEV, R2.DEV, na.rm=T)) / 100, ceiling(100 * max(est.R2.DEV, R2.DEV, na.rm=T)) / 100),
     main = "N=250")
grid()
abline(a = 0, b = 1, col = 'red')

plot(R2.DEV[,3], est.R2.DEV[,3], 
     pch  = 20, 
     xlab = expression(paste("model ", italic(R)) ["cs"]^"2"), 
     ylab = expression(paste("approximated ", italic(R)) ["cs"]^"2"), 
     ylim = c(floor(100 * min(est.R2.DEV, R2.DEV, na.rm=T)) / 100, ceiling(100 * max(est.R2.DEV, R2.DEV, na.rm=T)) / 100), 
     xlim = c(floor(100 * min(est.R2.DEV, R2.DEV, na.rm=T)) / 100, ceiling(100 * max(est.R2.DEV, R2.DEV, na.rm=T)) / 100),
     main = "N=500")
grid()
abline(a = 0, b = 1, col = 'red')

plot(R2.DEV[,4], est.R2.DEV[,4], 
     pch  = 20, 
     xlab = expression(paste("model ", italic(R)) ["cs"]^"2"), 
     ylab = expression(paste("approximated ", italic(R)) ["cs"]^"2"), 
     ylim = c(floor(100 * min(est.R2.DEV, R2.DEV, na.rm=T)) / 100, ceiling(100 * max(est.R2.DEV, R2.DEV, na.rm=T)) / 100), 
     xlim = c(floor(100 * min(est.R2.DEV, R2.DEV, na.rm=T)) / 100, ceiling(100 * max(est.R2.DEV, R2.DEV, na.rm=T)) / 100),
     main = "N=1000")
grid()
abline(a = 0, b = 1, col = 'red')

