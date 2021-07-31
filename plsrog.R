plsrog <- function(X,class, kappa){
  
  # data matrix
  X <- as.matrix(X)
  X <- matrix(as.numeric(X),nrow=nrow(X)) # metabolites*samples
  
  # response variable
  Y0 <- factor(class)
  Y <- model.matrix(~ Y0 + 0)
  
  # penalized matrix
  P <- NULL
  p <- colSums(Y)
  for(i in 1:ncol(Y)){
    P <- cbind(P,Y[,i]/p[i])
  }
  P <- t(P)
  
  # differential matrix
  g <- ncol(Y)
  D <- diff(diag(1,g))
    
  # autoscaling
  X <- scale(X)
  Y <- scale(Y,scale=FALSE)
  
  # sample size-1
  N <- nrow(X)-1
  
  # smoothing parameter
  C <- kappa*t(Y)%*%t(P)%*%t(D)%*%D%*%P%*%Y+(1-kappa)*diag(1,g)
  
  # cholesky decomposition
  Rx <- chol(solve(C))
  Ry <- chol(C)
  
  # singular value decomposition
  USVx <- svd(Rx%*%t(Y)%*%X/N)
  USVy <- svd(t(X)%*%Y%*%solve(Ry)/N)
  
  # weght vector (factor loading)
  Wx <- USVx$v
  Wy <- solve(Ry)%*%USVy$v
  
  # PLS score
  T <- X%*%Wx
  S <- Y%*%Wy
    
  # correlation coefficient
  R <- NULL
  for(i in 1:(ncol(Y))){
      lambdax <- cov(T[,i],S[,i])
      r <- (sqrt(N)*lambdax*Wx[,i])/as.numeric(sqrt(t(Wy[,i])%*%t(Y)%*%Y%*%Wy[,i]))
      r[is.nan(r)] <- 0
      R <- cbind(R,r)
  }
  
  # statistical test
  P <- NULL
  for(i in 1:(ncol(Y))){
    p <- 2*pt(abs(R[,i])*sqrt(nrow(X)-2)/sqrt(1-R[,i]^2), nrow(X)-2, lower.tail=FALSE)
    P <- cbind(P,p)
  }
  
  # Q value
  Q <- NULL
  for(i in 1:(ncol(Y))){
    q <- p.adjust(P[,i], method = "BH")
    # q <- -log10(q)
    Q <- cbind(Q,q)
  }
  
  all <- list(T,S,Wx,Wy,R,P,Q)
}

pls <- function(X,class){
  
  # data matrix
  X <- as.matrix(X)
  X <- matrix(as.numeric(X),nrow=nrow(X)) # metabolites*samples

  # response variable
  Y0 <- factor(class)
  Y <- model.matrix(~ Y0 + 0)
  
  # penalized matrix
  P <- NULL
  p <- colSums(Y)
  for(i in 1:ncol(Y)){
    P <- cbind(P,Y[,i]/p[i])
  }
  P <- t(P)
  
  # autoscaling
  X <- scale(X)
  Y <- scale(Y,scale=FALSE)
  
  # ----------------
  #   ordinary PLS
  # ----------------
  
  # (sample size)-1
  N <- nrow(X)-1
  
  # singular value decomposition
  USVx <- svd(t(Y)%*%X/N)
  USVy <- svd(t(X)%*%Y/N)
  
  # weight vector matrix
  Wx <- USVx$v
  Wy <- USVy$v
  
  # score matrix
  T <- X%*%Wx
  S <- Y%*%Wy
  
  # correlation coefficient
  R <- NULL
  for(i in 1:(ncol(Y))){
    lambdax <- cov(T[,i],S[,i])
    options(warn=-1)
    r <- (sqrt(N)*lambdax*Wx[,i])/(sqrt(t(Wy[,i])%*%t(Y)%*%Y%*%Wy[,i]))
    r[is.nan(r)] <- 0
    R <- cbind(R,r)
  }
  
  # statistical test
  P <- NULL
  for(i in 1:(ncol(Y))){
    p <- 2*pt(abs(R[,i])*sqrt(nrow(X)-2)/sqrt(1-R[,i]^2), nrow(X)-2, lower.tail=FALSE)
    P <- cbind(P,p)
  }
  
  # Q value
  Q <- NULL
  for(i in 1:(ncol(Y))){
    q <- p.adjust(P[,i], method = "BH")
    # q <- -log10(q)
    Q <- cbind(Q,q)
  }
  
  all <- list(T,S,Wx,Wy,R,P,Q)
}