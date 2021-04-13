dglm_laplace <- function(y, X, b0, Sig0, G, W, 
                         smoother = T, model){
  
  n <- length(y)
  b <- matrix(NA, nrow = n, ncol = length(b0))
  Sig <- array(NA, c(length(b0), length(b0), n))
  p <- rep(NA, n)
  
  # initialize
  b[1, ] <- b0
  Sig[, , 1] <- Sig0
  
  if(model == 'binomial'){
    p[1] <- 1/ (1 + exp(-t(X[1, ]) %*% b[1, ]))
  }else if(model == 'poisson'){
    p[1] <- exp(t(X[1, ]) %*% b[1, ])
  }else{
    cat('unsupported distribution')
    return()
  }
  
  bpred <- b
  Sigpred <- Sig
  
  # forward: filtering
  for(i in 2:n){
    bpred[i, ] <- G %*% b[i-1, ]
    Sigpred[, , i] <- G %*% Sig[, , i-1] %*% t(G) + W
    
    if(model == 'binomial'){
      
      p[i] <- 1/ (1 + exp(-t(X[i, ]) %*% bpred[i, ]))
      Sigpostinv <- solve(Sigpred[, , i]) +
        p[i]*(1-p[i])*X[i, ]%*%t(X[i, ])
      Sig[, , i] <- solve(Sigpostinv)
      b[i, ] <- bpred[i, ] + (y[i] - p[i])* Sig[, , i] %*% X[i, ]
      
    }else if(model == 'poisson'){
      p[i] <- exp(t(X[i, ]) %*% bpred[i, ])
      Sigpostinv <- solve(Sigpred[, , i]) +
        p[i]*X[i, ]%*%t(X[i, ])
      Sig[, , i] <- solve(Sigpostinv)
      b[i, ] <- bpred[i, ] + (y[i] - p[i])* Sig[, , i] %*% X[i, ]
    }
    
  }
  
  # backward: smoothing
  if(smoother){
    for(j in (n-2):1){
      C <- Sig[, , j] %*% G %*% solve(Sigpred[, , j+1])
      b[j, ] <- b[j, ] + C %*% (b[j+1, ] - bpred[j+1, ])
      Sig[, , j] <- Sig[, , j] + C %*% (Sig[, , j+1] - Sigpred[, , j+1]) %*% t(C)
      
    }
  }
  return(list(b= b, Sig= Sig, p = p))
  
}


dglm_conBayes <- function(y, X, b0, Sig0, G, W,
                          smoother = T, model){
  
  n <- length(y)
  b <- matrix(NA, nrow = n, ncol = length(b0))
  Sig <- array(NA, c(length(b0), length(b0), n))
  
  # initialize
  b[1, ] <- b0
  Sig[, , 1] <- Sig0
  
  bpred <- b
  Sigpred <- Sig
  
  for(i in 2:n){
    bpred[i, ] <- G %*% b[i-1, ]
    Sigpred[, , i] <- G %*% Sig[, , i-1] %*% t(G) + W
    st <- Sigpred[, , i] %*% X[i, ]
    
    ft <- t(X[i, ]) %*% bpred[i, ]
    qt <- t(X[i, ]) %*% st
    
    if(model == 'binomial'){
      at <- (1/qt)*(1 + exp(ft))
      bt <- (1/qt)*(1 + exp(-ft))
      gt <- digamma(at + y[i]) - digamma(bt + 1 - y[i])
      pt <- trigamma(at + y[i]) + trigamma(bt + 1 - y[i])
      
    }else if(model == 'poisson'){
      at <- 1/qt
      bt <- (1/qt)*exp(-ft)
      gt <- digamma(at + y[i]) - log(bt + 1)
      pt <- trigamma(at + y[i])
      
    }
    
    b[i, ] <- bpred[i, ] + st *c((gt - ft)/qt)
    Sig[, , i] <- Sigpred[, , i] - st %*% t(st) * c((1 - pt/qt)/qt)
    
  }
  
  # backward: smoothing
  if(smoother){
    for(j in (n-2):1){
      C <- Sig[, , j] %*% G %*% solve(Sigpred[, , j+1])
      b[j, ] <- b[j, ] + C %*% (b[j+1, ] - bpred[j+1, ])
      Sig[, , j] <- Sig[, , j] + C %*% (Sig[, , j+1] - Sigpred[, , j+1]) %*% t(C)
      
    }
  }
  return(list(b= b, Sig= Sig))
  
}

#############################################
############ Poisson Example ################
#############################################

######## simulation
set.seed(2)
n <- 1000
period <- n/3
X <- cbind(rnorm(n, 1, 1), rnorm(n, 1, 1))
beta <- cbind(
  c(rep(0, round(n/2)), rep(2, n - round(n/2))),
  1 +sin((2*pi/period)*(1:n))
)

eta <- rowSums(X*beta)
p <- exp(eta)
y <- rpois(n, p)

b0 <- solve(t(X[1:2, ]), rep(log(mean(y[1:100])), 2))
G <- diag(1, 2)
Sig0 <- diag(1, 2)*1e-1
W <- diag(1, 2)*1e-3

######## fitting
## Laplace approximation
# filtering
res_LF <- dglm_laplace(y, X, b0, Sig0, G, W, 
                       smoother = F, model = 'poisson')

# smoothing
res_LS <- dglm_laplace(y, X, b0, Sig0, G, W, 
                       smoother = T, model = 'poisson')

## conjugate prior
# filtering
res_CF <- dglm_conBayes(y, X, b0, Sig0, G, W,
                        smoother = F, model = 'poisson')

# smoothing
res_CS <- dglm_conBayes(y, X, b0, Sig0, G, W,
                        smoother = T, model = 'poisson')

######## deviance
deviance <- function(res){
  etaPred <- rowSums(X*res$b)
  pPred <- 1/ (1 + exp(-etaPred))
  return(2*sum(ifelse(y == 0, 0, y*log(y/pPred)) - (y - pPred)))
}

deviance(res_LF)
deviance(res_LS)
deviance(res_CF)
deviance(res_CS)

######## plot
#### Laplace approximation
plot(beta[, 1], type = 'l', lwd = 3,
     ylim = c(min(beta[, 1])-.6, max(beta[, 1])+.3),
     xlab = 't', ylab = 'first state vector')
lines(res_LF$b[, 1], col = 'red', lwd = 3)
lines(res_LF$b[, 1] + res_LF$Sig[1, 1, ], col = 'red', lty = 2, lwd = 3)
lines(res_LF$b[, 1] - res_LF$Sig[1, 1, ], col = 'red', lty = 2, lwd = 3)
lines(res_LS$b[, 1], col = 'blue', lwd = 3)
lines(res_LS$b[, 1] + res_LF$Sig[1, 1, ], col = 'blue', lty = 2, lwd = 3)
lines(res_LS$b[, 1] - res_LF$Sig[1, 1, ], col = 'blue', lty = 2, lwd = 3)
legend('topleft', legend = c('true', 'filtering', 'smoothing'),
       lwd = 3, col = c('black', 'red', 'blue'))

plot(beta[, 2], type = 'l', lwd = 3,
     ylim = c(min(beta[, 2])-.3, max(beta[, 2])+.3),
     xlab = 't', ylab = 'second state vector')
lines(res_LF$b[, 2], col = 'red', lwd = 3)
lines(res_LF$b[, 2] + res_LF$Sig[2, 2, ], col = 'red', lty = 2, lwd = 3)
lines(res_LF$b[, 2] - res_LF$Sig[2, 2, ], col = 'red', lty = 2, lwd = 3)
lines(res_LS$b[, 2], col = 'blue', lwd = 3)
lines(res_LS$b[, 2] + res_LF$Sig[2, 2, ], col = 'blue', lty = 2, lwd = 3)
lines(res_LS$b[, 2] - res_LF$Sig[2, 2, ], col = 'blue', lty = 2, lwd = 3)
legend('topleft', legend = c('true', 'filtering', 'smoothing'),
       lwd = 3, col = c('black', 'red', 'blue'))

#### conjugate prior
plot(beta[, 1], type = 'l', lwd = 3,
     ylim = c(min(beta[, 1])-.6, max(beta[, 1])+.3),
     xlab = 't', ylab = 'first state vector')
lines(res_CF$b[, 1], col = 'red', lwd = 3)
lines(res_CF$b[, 1] + res_CF$Sig[1, 1, ], col = 'red', lty = 2, lwd = 3)
lines(res_CF$b[, 1] - res_CF$Sig[1, 1, ], col = 'red', lty = 2, lwd = 3)
lines(res_CS$b[, 1], col = 'blue', lwd = 3)
lines(res_CS$b[, 1] + res_CF$Sig[1, 1, ], col = 'blue', lty = 2, lwd = 3)
lines(res_CS$b[, 1] - res_CF$Sig[1, 1, ], col = 'blue', lty = 2, lwd = 3)
legend('topleft', legend = c('true', 'filtering', 'smoothing'),
       lwd = 3, col = c('black', 'red', 'blue'))

plot(beta[, 2], type = 'l', lwd = 3,
     ylim = c(min(beta[, 2])-.3, max(beta[, 2])+.3),
     xlab = 't', ylab = 'second state vector')
lines(res_CF$b[, 2], col = 'red', lwd = 3)
lines(res_CF$b[, 2] + res_CF$Sig[2, 2, ], col = 'red', lty = 2, lwd = 3)
lines(res_CF$b[, 2] - res_CF$Sig[2, 2, ], col = 'red', lty = 2, lwd = 3)
lines(res_CS$b[, 2], col = 'blue', lwd = 3)
lines(res_CS$b[, 2] + res_CF$Sig[2, 2, ], col = 'blue', lty = 2, lwd = 3)
lines(res_CS$b[, 2] - res_CF$Sig[2, 2, ], col = 'blue', lty = 2, lwd = 3)
legend('topleft', legend = c('true', 'filtering', 'smoothing'),
       lwd = 3, col = c('black', 'red', 'blue'))


#############################################
############ Binomial Example ###############
#############################################

######## simulation
set.seed(1)
n <- 1000
X <- cbind(rnorm(n, 1, 1), rnorm(n, 1, 1))
beta <- cbind(
  seq(0, 2, length.out = n),
  seq(2, 0, length.out = n) #
)

eta <- rowSums(X*beta)
p <- 1/ (1 + exp(-eta))
y <- rbinom(n, 1, p)

b0 <- beta[1, ]
G <- diag(1, 2)
Sig0 <- diag(1, 2)*1e-1
W <- diag(1, 2)*1e-3

######## fitting
# filtering
res_LF <- dglm_laplace(y, X, b0, Sig0, G, W, 
                       smoother = F, model = 'binomial')

# smoothing
res_LS <- dglm_laplace(y, X, b0, Sig0, G, W, 
                       smoother = T, model = 'binomial')


######## deviance
deviance(res_LF)
deviance(res_LS)


######## plot
plot(beta[, 1], type = 'l', lwd = 3,
     ylim = c(min(beta[, 1])-.3, max(beta[, 1])+.3),
     xlab = 't', ylab = 'first state vector')
lines(res_LF$b[, 1], col = 'red', lwd = 3)
lines(res_LF$b[, 1] + res_LF$Sig[1, 1, ], col = 'red', lty = 2, lwd = 3)
lines(res_LF$b[, 1] - res_LF$Sig[1, 1, ], col = 'red', lty = 2, lwd = 3)
lines(res_LS$b[, 1], col = 'blue', lwd = 3)
lines(res_LS$b[, 1] + res_LF$Sig[1, 1, ], col = 'blue', lty = 2, lwd = 3)
lines(res_LS$b[, 1] - res_LF$Sig[1, 1, ], col = 'blue', lty = 2, lwd = 3)
legend('topleft', legend = c('true', 'filtering', 'smoothing'),
       lwd = 3, col = c('black', 'red', 'blue'))

plot(beta[, 2], type = 'l', lwd = 3,
     ylim = c(min(beta[, 2])-.3, max(beta[, 2])+.3),
     xlab = 't', ylab = 'second state vector')
lines(res_LF$b[, 2], col = 'red', lwd = 3)
lines(res_LF$b[, 2] + res_LF$Sig[2, 2, ], col = 'red', lty = 2, lwd = 3)
lines(res_LF$b[, 2] - res_LF$Sig[2, 2, ], col = 'red', lty = 2, lwd = 3)
lines(res_LS$b[, 2], col = 'blue', lwd = 3)
lines(res_LS$b[, 2] + res_LF$Sig[2, 2, ], col = 'blue', lty = 2, lwd = 3)
lines(res_LS$b[, 2] - res_LF$Sig[2, 2, ], col = 'blue', lty = 2, lwd = 3)
legend('topleft', legend = c('true', 'filtering', 'smoothing'),
       lwd = 3, col = c('black', 'red', 'blue'))

#############################################










