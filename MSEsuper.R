library(HHG)
library(minerva) 
library(energy)
library(GeneralizedHyperbolic)
library(ggplot2)

# 普通光滑
MSEn <- function(model, n, sd0, a) { 
  
  # 设置噪声
  if (model <= 2) {
    noise <- rnorm(n)
  } else {
    noise <- rep(1, n)
  }
  
  repeat {
  # 生成数据根据不同模型
  if (model == 1) {  # 线性模型
    xs <- rnorm(n)
    ys <- a * xs + noise
  } else if (model == 2) {  # 二次模型
    xs <- rnorm(n)
    ys <- a * xs^2 + noise
  } 
  
  # 生成可观测数据
  us <- matrix(rnorm(n, mean = 0, sd = sd0), ncol = 1) 
  ws <- xs + us  # 可观测数据
  
  X <- ws  # 输入X替换为可观测的ws
  Y <- ys
  
  # 计算统计量
  source("select_h.R")
  h1_values <- seq(0.1, 1, by = 0.01) 
  h2_values <- seq(0.5, 3, by = 0.03) 
  #设置窗宽
  HH <- besth1(sd0, h1_values, h2_values, X, Y, a, noise)
  h1 <- HH$best_h1 
  h2 <- HH$best_h2
  h3 <- h2
  # 初始化变量
  f0 <- rep(n)
  f1 <- rep(n)
  f2 <- rep(n)
  f3 <- rep(n)
  # 计算f0, f1, f2, f3
  for (i in 1:n) { 
    f0[i] <- fdc(Y[i], Y[-i], h3)
    f1[i] <- fDC(sd0, X[i], Y[i], X[-i], Y[-i], h1, h2) / f0[i]
    f2[i] <- fNW(X[i], Y[i], X[-i], Y[-i], h1, h2) / f0[i]
    f3[i] <- fLL(X[i], Y[i], X[-i], Y[-i], h1, h2) / f0[i]
  }
  
  # 定义一个函数来处理列
  remove_invalid_values <- function(col) {
    valid_values <- col[is.finite(col) & col > 0]
    return(valid_values)
  }
  
  valid_pp1 <- remove_invalid_values(f1)
  f01 <- 1 / n * sum(log(valid_pp1))
  
  valid_pp2 <- remove_invalid_values(f2)
  f02 <- 1 / n * sum(log(valid_pp2))
  
  valid_pp3 <- remove_invalid_values(f3)
  f03 <- 1 / n * sum(log(valid_pp3))
  
  if (abs(f01) <= 1) {
    break
  }
}

  # 返回计算得到的统计量
  return(list(f01 = f01, f02 = f02, f03 = f03))
  #print(return(list(f01 = f01, f02 = f02, f03 = f03)))
}


# 定义进行N次循环的函数
calculate_MSE <- function(N, model, n, sd0, a) {
  # 初始化变量用于存储每次循环的结果
  fDC_diffs <- numeric(N)
  fNW_diffs <- numeric(N)
  fLL_diffs <- numeric(N)
  
  for (i in 1:N) {
    result <- MSEn(model, n, sd0, a)
    #print(result)
    fDC_diffs[i] <- result$f01 
    fNW_diffs[i] <- result$f02
    fLL_diffs[i] <- result$f03
  }
  
  # 计算均方误差
  MSE_fDC <- sum((fDC_diffs - mean(fDC_diffs))^2) / N
  MSE_fNW <- sum((fNW_diffs- mean(fNW_diffs))^2) / N
  MSE_fLL <- sum((fLL_diffs- mean(fLL_diffs))^2) / N
  
  return(list(MSE_fDC = MSE_fDC, MSE_fNW = MSE_fNW, MSE_fLL = MSE_fLL))
}

# 设置参数
N <- 50
model <- 2
#ns <- c(50, 100, 200, 500, 1000)
n <-100
sd0 <- sqrt(0.25)
a <- 0
print(calculate_MSE(N, model, n, sd0, a))

