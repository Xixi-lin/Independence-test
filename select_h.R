#model1：Y=aX+b (b=noise)
# 定义概率密度函数f(y) 
f_y1 <- function(y, a, b) {
  return (1 / sqrt(2 * pi * (a^2 + b^2))) * exp(-y^2 / (2 * (a^2 + b^2)))
}

# 定义条件概率密度函数f(y|x)
fy_given_x1 <- function(y, x, a, b) {
  return (1 / sqrt(2 * pi * b^2)) * exp(-((y - a * x)^2) / (2 * b^2))
}

#model2
# 定义概率密度函数f(y)
f_y2 <- function(y, a, b) {
  integrand <- function(x) {
    exp_val <- exp(-x^2 / 2)
    e_val <- exp(-((y - a * x^2) / b)^2 / 2) / b
    return(exp_val * e_val * abs(x))
  }
  integral_result <- integrate(integrand, -Inf, Inf)$value
  return(integral_result)
}

# 定义条件概率密度函数f(y|x)
fy_given_x2 <- function(y, x, a, b) {
  return (exp(-((y - a * x^2) / b)^2 / 2) / (sqrt(2 * pi) * b))
}

#model3
# 定义概率密度函数f(y) 

# 定义条件概率密度函数f(y|x)


#super smooth ker
Kn1=function(x,sigma_0,h_n){
  f <- function(t) {
    return(1/pi*cos(t*x)*(1-t^2)^3*exp((sigma_0^2*t^2)/(2*h_n^2))) #u为正态分布
  }
  result <- integrate(f, lower = 0, upper = 1)
  integral_value <- result$value
  return(integral_value)
}

#ordinary smooth ker
Kn2=function(x,sigma_0,h_n){
  return(1/sqrt(2*pi)*exp(-1/2*x^2)*(1-sigma_0^2/(2*h_n^2)*(x^2-1))) #u为双指数分布
}

#classic ker
ker <- function(y){
  k = 3/4*(1-y^2)*(abs(y)<=1)
  # k=dnorm(y)
  return(k) 
}

#DC method:f(yi|xi)
fDC <- function(sd0, w, y, W, Y, h1, h2){ #W,Y要移除第i个(w,y)
    m <- length(W) #n-1
    a2=matrix(0,m,1)
    for (j in 1:m){
      #a2[j]=Kn1((w-W[j])/h1,sd0,h1) #super smooth正态分布
      a2[j]=Kn2((w-W[j])/h1,sd0,h1) #ordinary smooth双指数分布
    }
    G2=sum(a2)
    K=ker((y-Y)/h2)/h2
    GK2=t(a2)%*%K
    fyx=GK2/G2
    return(fyx)
}

#DC method:f(yi),与NW, LL相同
fdc <- function(y, Y, h3){ #Y要移除第i个(y)
  n <- length(Y)
  fy=sum(ker((y-Y)/h3))/((n-1)*h3)
  return(fy) 
}

#NW method:f(yi|xi)
fNW <- function(w, y, W, Y, h1, h2){ #W,Y要移除第i个(w,y)
  a1=ker((w-W)/h1)
  G1=sum(a1)
  K=ker((y-Y)/h2)/h2
  GK1=t(a1)%*%K
  fyx=GK1/G1
  return(fyx)
}

#LL method:f(yi|xi)
fLL <- function(w, y, W, Y, h1, h2){ #W,Y要移除第i个(w,y)
  m0=sum(ker((w-W)/h1))
  m1=sum(ker((w-W)/h1)*(W-w))
  m2=sum(ker((w-W)/h1)*(W-w)^2)
  a3=m0*m2-m1^2
  G3=sum(a3)
  K=ker((y-Y)/h2)/h2
  GK3=t(ker((w-W)/h1)*(m2-(W-w)*m1))%*%K
  fyx=GK3/G3
  return(fyx)
}

#model1
Besth1 <- function(h3_values, Y, a, b) {
  num_h3 <- length(h3_values)
  n <- length(Y)
  r_values <- numeric(num_h3)  # Initialize a numeric vector to store r values
  for (i in 1:num_h3) {
    h3 <- h3_values[i]
    fy <- numeric(n)
    for (xx in 1:n) {
      fy[xx]<-fdc(Y[xx],Y[-xx],h3)
    }
    fy0 = f_y1(Y, a, b)
    r <- mean((fy - fy0)^2)
    r_values[i] <- r  # Store the r value in the vector
  }
  min_r <- min(r_values)  # Find the minimum r value
  best_h3 <- h3_values[which.min(r_values)]
  return(best_h3)
}

besth1 <- function(sd0, h1_values, h2_values, W, Y, a, b) {
  num_h1 <- length(h1_values)
  num_h2 <- length(h2_values)
  n <- length(Y)
  r2_values <- matrix(0, nrow = num_h1, ncol = num_h2)  # Initialize matrix to store r2 values
  for (i in 1:num_h1) {
    h1 <- h1_values[i]
    for (j in 1:num_h2) {
      h2 <- h2_values[j]
      fyx <- numeric(n)
      for (k in 1:n) {
        fyx[k] <- fDC(sd0, W[k], Y[k], W[-k], Y[-k], h1, h2)
      }
      fyx0 <- fy_given_x1(Y, W, a, b)
      r2 <- mean((fyx-fyx0)^2)
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

#model2
Besth2 <- function(h3_values, Y, a, b) {
  num_h3 <- length(h3_values)
  n <- length(Y)
  r_values <- numeric(num_h3)  # Initialize a numeric vector to store r values
  for (i in 1:num_h3) {
    h3 <- h3_values[i]
    fy <- numeric(n)
    for (xx in 1:n) {
      fy[xx]<-fdc(Y[xx],Y[-xx],h3)
    }
    fy0 <- f_y2(Y, a, b) 
    r <- mean((fy - fy0)^2)
    r_values[i] <- r  # Store the r value in the vector
  }
  min_r <- min(r_values)  # Find the minimum r value
  best_h3 <- h3_values[which.min(r_values)]
  return(best_h3)
}

besth2 <- function(sd0, h1_values, h2_values, W, Y, a, b) {
  num_h1 <- length(h1_values)
  num_h2 <- length(h2_values)
  n <- length(Y)
  r2_values <- matrix(0, nrow = num_h1, ncol = num_h2)  # Initialize matrix to store r2 values
  for (i in 1:num_h1) {
    h1 <- h1_values[i]
    for (j in 1:num_h2) {
      h2 <- h2_values[j]
      fyx <- numeric(n)
      for (k in 1:n) {
        fyx[k] <- fDC(sd0, W[k], Y[k], W[-k], Y[-k], h1, h2)
      }
      fyx0 <- fy_given_x2(Y, W, a, b)
      #fyx0 <- dnorm(Y, mean = muY, sd = sqrt(SigY))
      r2 <- mean((fyx-fyx0)^2)
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

#model3, model4
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


