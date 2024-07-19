library(HHG)
library(minerva) 
library(energy)
library(GeneralizedHyperbolic)
library(ggplot2)

# 普通光滑
MSEd <- function(model, n, sd0, a) { 
  
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
  mu <- 0
  beta <- sd0 / sqrt(2)
  alpha <- beta
  us <- rskewlap(n, mu, alpha, beta, param = c(mu, alpha, beta))
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
  #adcor <- dcor(X, Y)
  #amic <- mine(X, Y)$MIC
  #Dx <- as.matrix(dist(X), diag = TRUE, upper = TRUE)
  #Dy <- as.matrix(dist(Y), diag = TRUE, upper = TRUE)
  # 计算 HHG 统计量
  #perm.pval.hhg.sc #Sum of chi-squared scores（卡方检验的总和）的置换检验的p值
  #perm.pval.hhg.sl #Sum of likelihood ratio scores（似然比检验的总和）的置换检验的p值
  #perm.pval.hhg.mc #Maximum over chi-squared scores（卡方检验的最大值）的置换检验的p值
  #perm.pval.hhg.ml #Maximum over likelihood ratio scores（似然比检验的最大值）的置换检验的p值
  #aHHG <- hhg.test(Dx, Dy)
  #aHHG <- aHHG$perm.pval.hhg.sc
 
  
  # 返回计算得到的统计量
  return(list(f01 = f01, f02 = f02, f03 = f03))
}

# 调用函数并提取结果
#result <- MSEd(1, 100, sqrt(0.1), 0)
#print(result)

# 定义进行N次循环的函数
calculate_MSE <- function(N, model, n, sd0, a) {
  # 初始化变量用于存储每次循环的结果
  #adcor_diffs <- numeric(N)
  #amic_diffs <- numeric(N)
  #aHHG_diffs <- numeric(N)
  fDC_diffs <- numeric(N)
  fNW_diffs <- numeric(N)
  fLL_diffs <- numeric(N)
  
  for (i in 1:N) {
    result <- MSEd(model, n, sd0, a)
    #adcor_diffs[i] <- (result$adcor - mean(result$adcor))^2
    #amic_diffs[i] <- (result$amic - mean(result$amic0))^2
    #aHHG_diffs[i] <- (result$aHHG - mean(result$aHHG0))^2
    fDC_diffs[i] <- result$f01 
    fNW_diffs[i] <- result$f02
    fLL_diffs[i] <- result$f03
  }
  
  # 计算均方误差
  #MSE_adcor <- sum(adcor_diffs) / N
  #MSE_amic <- sum(amic_diffs) / N
  #MSE_aHHG <- sum(aHHG_diffs) / N
  MSE_fDC <- sum((fDC_diffs - mean(fDC_diffs))^2) / N
  MSE_fNW <- sum((fNW_diffs- mean(fNW_diffs))^2) / N
  MSE_fLL <- sum((fLL_diffs- mean(fLL_diffs))^2) / N
  
  return(list(MSE_fDC = MSE_fDC, MSE_fNW = MSE_fNW, MSE_fLL = MSE_fLL))
}

# 设置参数
N <- 30
model <- 1
#ns <- c(50, 100, 200, 500, 1000)
n <-100
sd0 <- sqrt(0.1)
a <- 0.6
print(calculate_MSE(N, model, n, sd0, a))

# 初始化数据框
#mse_results_df <- data.frame(n = integer(),
#                             method = character(),
#                             MSE = numeric())

# 进行计算并存储结果
#for (n in ns) {
#  mse_results <- calculate_MSE(N, model, n, sd0, a)
#  mse_results_df <- rbind(mse_results_df, data.frame(n = n, method = "adcor", MSE = mse_results$MSE_adcor))
#  mse_results_df <- rbind(mse_results_df, data.frame(n = n, method = "amic", MSE = mse_results$MSE_amic))
#  mse_results_df <- rbind(mse_results_df, data.frame(n = n, method = "aHHG", MSE = mse_results$MSE_aHHG))
#  mse_results_df <- rbind(mse_results_df, data.frame(n = n, method = "fDC", MSE = mse_results$MSE_fDC))
#  mse_results_df <- rbind(mse_results_df, data.frame(n = n, method = "fNW", MSE = mse_results$MSE_fNW))
#  mse_results_df <- rbind(mse_results_df, data.frame(n = n, method = "fLL", MSE = mse_results$MSE_fLL))
#}

#print(mse_results_df)


# 绘制图形
#p <- ggplot(mse_results_df, aes(x = n, y = MSE, color = method, linetype = method)) +
#  geom_smooth(se = FALSE, size = 0.7, linetype = "solid") + # 添加平滑曲线
#  scale_color_manual(values = c("adcor" = "black", "amic" = "brown", 
#                                "aHHG" = "purple", "fDC" = "red", 
#                                "fNW" = "blue", "fLL" = "green"),
#                     breaks = c("fDC", "fNW", "fLL", "amic", "adcor", "aHHG"),
#                     labels = c("DCT", "NWT", "LLT", "MIC", "dCor", "HHG")) +
#  scale_linetype_manual(values = c("adcor" = "dotdash", "amic" = "longdash", 
#                                   "aHHG" = "twodash", "fDC" = "solid", 
#                                   "fNW" = "dashed", "fLL" = "dotted"),
#                        breaks = c("fDC", "fNW", "fLL", "amic", "adcor", "aHHG"),
#                        labels = c("DCT", "NWT", "LLT", "MIC", "dCor", "HHG")) +
#  labs(x = "Sample Size", y = "MSE", color = "", linetype = "") +
#  theme_minimal() +
#  theme(
#    panel.border = element_rect(color = "black", fill = NA, size = 1), # 添加边框
#    panel.grid.major = element_blank(), # 去掉主要网格线
#    panel.grid.minor = element_blank(),  # 去掉次要网格线
#    legend.position = c(0.1, 0.80), # 设置图例位置为右上角
#    legend.background = element_rect(fill = "white", color = "black"),
#    legend.key.width = unit(0.8, "cm"),  # 设置图例框的宽度为1厘米
#    legend.key.height = unit(0.5, "cm"),  # 设置图例框的高度为0.5厘米
#    legend.margin = margin(-15.0, 0.8, 0.8, 1.0)  # 设置图例的边距，上、右、下、左
#  )

#print(p)




