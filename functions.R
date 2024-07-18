Kn1=function(x,sigma_0,h_n){
  f <- function(t) {
    return(1/pi*cos(t*x)*(1-t^2)^3*exp((sigma_0^2*t^2)/(2*h_n^2))) #u为正态分布
  }
  result <- integrate(f, lower = 0, upper = 1)
  integral_value <- result$value
  return(integral_value)
}

Kn2=function(x,sigma_0,h_n){
    return(1/sqrt(2*pi)*exp(-1/2*x^2)*(1-sigma_0^2/(2*h_n^2)*(x^2-1))) #u为双指数分布
}

ker <- function(y){
  k = 3/4*(1-y^2)*(abs(y)<=1)
  # k=dnorm(y)
  return(k) 
}

#example1\change-e
besth1 <- function(x0, muX,muY,SigX,SigY, sd0, rho, h1_values, h2_values, Y, W, y) {
  num_h1 <- length(h1_values)
  num_h2 <- length(h2_values)
  M <- length(y)
  n <- length(Y)
  
  r2_values <- matrix(0, nrow = num_h1, ncol = num_h2)  # Initialize matrix to store r2 values
  
  for (i in 1:num_h1) {
    h1 <- h1_values[i]
    for (j in 1:num_h2) {
      h2 <- h2_values[j]
      
      GK2 <- matrix(0,M,1)
      f2 <- matrix(0,M,1)
      
      a2=matrix(0,n,1)
      for (k in 1:n) {
        a2[k] <- Kn1((x0 - W[k]) / h1, sd0, h1) #1
        #a2[k] <- Kn2((x0 - W[k]) / h1, sd0, h1) #change-e
      }
      G2=sum(a2)
      
      for (xx in 1:M) {
        K2=ker((Y-y[xx])/h2)/h2
        GK2[xx]=t(a2)%*%matrix(K2,n,1)
        f2[xx] <- GK2[xx] / G2
      }
      
      f0 <- 1/sqrt(2*pi*(1-rho^2)*SigY)*exp(-(rho*(x0-muX)/sqrt(SigX)-(y-muY)/sqrt(SigY))^2/(2*(1-rho^2)))
      r2 <- sqrt(mean((f2-f0)^2))
      
      r2_values[i, j] <- r2  # Store the r2 value in the matrix
    }
  }
  
  min_r2 <- min(r2_values)  # Find the minimum r2 value
  ind <- which.min(r2_values)  # Find the indices of the minimum r2 value
  loc <- arrayInd(ind, dim(r2_values))
  best_h1 <- h1_values[loc[1]]
  best_h2 <- h2_values[loc[2]]
  
  return(list(best_h1 = best_h1, best_h2 = best_h2))
}



#example2.1 2.2
besth2 <- function(x0, sd0, h1_values, h2_values, Y, W, y, c, d) {
  num_h1 <- length(h1_values)
  num_h2 <- length(h2_values)
  M <- length(y)
  n <- length(Y)
  
  r2_values <- matrix(0, nrow = num_h1, ncol = num_h2)  # Initialize matrix to store r2 values
  
  for (i in 1:num_h1) {
    h1 <- h1_values[i]
    for (j in 1:num_h2) {
      h2 <- h2_values[j]
      
      GK2 <- matrix(0,M,1)
      f2 <- matrix(0,M,1)
      
      a2=matrix(0,n,1)
      for (k in 1:n) {
        a2[k] <- Kn1((x0-W[k])/h1, sd0, h1) #2.1
        #a2[k] <- Kn2((x0-W[k])/h1, sd0, h1) #2.2
      }
      G2=sum(a2)
      
      for (xx in 1:M) {
        K2=ker((Y-y[xx])/h2)/h2
        GK2[xx]=t(a2)%*%matrix(K2,n,1)
        f2[xx] <- GK2[xx] / G2
      }
      
      f0 <- 1/(c*sqrt(2*pi))*exp(-(y-d*exp(x0)-sin(pi*x0))^2/(2*c^2)) #s[]系数为1
      r2 <- sqrt(mean((f2-f0)^2))
      
      r2_values[i, j] <- r2  # Store the r2 value in the matrix
    }
  }
  
  min_r2 <- min(r2_values)  # Find the minimum r2 value
  ind <- which.min(r2_values)  # Find the indices of the minimum r2 value
  loc <- arrayInd(ind, dim(r2_values))
  best_h1 <- h1_values[loc[1]]
  best_h2 <- h2_values[loc[2]]
  
  return(list(best_h1 = best_h1, best_h2 = best_h2))
}


#example3.1 3.2
besth3 <- function(x0, sd0, h1_values, h2_values, Y, W, y) {
  num_h1 <- length(h1_values)
  num_h2 <- length(h2_values)
  M <- length(y)
  n <- length(Y)
  
  r2_values <- matrix(0, nrow = num_h1, ncol = num_h2)  # Initialize matrix to store r2 values
  
  for (i in 1:num_h1) {
    h1 <- h1_values[i]
    for (j in 1:num_h2) {
      h2 <- h2_values[j]
      
      GK2 <- matrix(0,M,1)
      f2 <- matrix(0,M,1)
      
      a2=matrix(0,n,1)
      for (k in 1:n) {
        #a2[k] <- Kn1((x0-W[k])/h1, sd0, h1) #3.1
        a2[k] <- Kn2((x0-W[k])/h1, sd0, h1) #3.2
      }
      G2=sum(a2)
      for (xx in 1:M) {
        K2=ker((Y-y[xx])/h2)/h2
        GK2[xx]=t(a2)%*%matrix(K2,n,1)
        f2[xx] <- GK2[xx]/G2
      }
      
      f0=1/(sqrt(2*pi)*0.5)*exp(-(y-x0)^2/(2*0.5^2))
      r2 <- sqrt(mean((f2-f0)^2))
      
      r2_values[i, j] <- r2  # Store the r2 value in the matrix
    }
  }
  
  min_r2 <- min(r2_values)  # Find the minimum r2 value
  ind <- which.min(r2_values)  # Find the indices of the minimum r2 value
  loc <- arrayInd(ind, dim(r2_values))
  best_h1 <- h1_values[loc[1]]
  best_h2 <- h2_values[loc[2]]
  
  return(list(best_h1 = best_h1, best_h2 = best_h2))
}



#太湖站:real data
#f(y|x)未知的h选择方式(#normal)
optimalh1 <- function(sd0, h1_values, h2_values, Y, W) {
  num_h1 <- length(h1_values)
  num_h2 <- length(h2_values)
  n <- length(Y)
  
  L_values <- matrix(0, nrow = num_h1, ncol = num_h2)  # Initialize matrix to store L values
  
  for (i in 1:num_h1) {
    h1 <- h1_values[i]
    
    for (j in 1:num_h2) {
      h2 <- h2_values[j]
      
      a4 <- matrix(0, n-1, n)
      K4 <- matrix(0, n-1, n)  # Adjust the dimensions of K2
      GK4 <- rep(0, n)      
      
      for (t in 1:n) {
        for (k in 1:(n-1)) {
          # if (k != t) {
          WW=W[-t]
          a4[k, t] <- Kn1((WW[k] - W[t])/h1, sd0, h1)
          # }
        }
        
        Y_diff <- Y[-t] - Y[t]
        K4[, t] <- ker(Y_diff / h2) / h2
        
        GK4[t] <- t(a4[ ,t]) %*% K4[ ,t] # Adjust the calculation for GK2
        GK4[t] <- (abs(GK4[t])/((n-1)*h1))^(1/n)
      }
      L <- prod(GK4)
      L_values[i, j] <- -L
    }
  }
  
  min_L <- min(L_values)
  ind <- which.min(L_values)
  loc <- arrayInd(ind, dim(L_values))
  optimal_h1 <- h1_values[loc[1]]
  optimal_h2 <- h2_values[loc[2]]
  
  return(list(optimal_h1 = optimal_h1, optimal_h2 = optimal_h2))
}

#f(y|x)未知的h选择方式(#double)
optimalh2 <- function(sd0, h1_values, h2_values, Y, W) {
  num_h1 <- length(h1_values)
  num_h2 <- length(h2_values)
  n <- length(Y)
  
  L_values <- matrix(0, nrow = num_h1, ncol = num_h2)  # Initialize matrix to store L values
  
  for (i in 1:num_h1) {
    h1 <- h1_values[i]
    
    for (j in 1:num_h2) {
      h2 <- h2_values[j]
      
      a4 <- matrix(0, n-1, n)
      K4 <- matrix(0, n-1, n)  # Adjust the dimensions of K2
      GK4 <- rep(0, n)      
      
      for (t in 1:n) {
        WW=W[-t]
        a4[ ,t] <- Kn2((WW - W[t])/h1, sd0, h1)
        Y_diff <- Y[-t] - Y[t]
        K4[ ,t] <- ker(Y_diff / h2) / h2
        
        GK4[t] <- t(a4[ ,t]) %*% K4[ ,t] # Adjust the calculation for GK2
        GK4[t] <- (abs(GK4[t])/((n-1)*h1))^(1/n)
      }
      
      L <- prod(GK4)
      L_values[i, j] <- -L
  
    }
  }
  
  min_L <- min(L_values)
  ind <- which.min(L_values)
  loc <- arrayInd(ind, dim(L_values))
  optimal_h1 <- h1_values[loc[1]]
  optimal_h2 <- h2_values[loc[2]]
  
  return(list(optimal_h1 = optimal_h1, optimal_h2 = optimal_h2))
}

