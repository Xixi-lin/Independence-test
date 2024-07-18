#从X,Y的关系入手生成数据
rm(list=ls())
time1<-Sys.time()
library(MASS)
library(ggplot2)
library(ggpubr)
library(latex2exp)
library(GeneralizedHyperbolic)
source("functions.R")

x0=1
sd0=sqrt(0.25) #sigma_0; b=sd0/sqrt(2)
c=0.6 #epsilon的系数
d=0.5
n<- 200 #样本数
y=seq(-0.5, 3, by = 0.1)
M=length(y) #y点个数
N=100 #循环次数
#h1_values <- seq(0.05, 0.6, by = 0.01)
#h2_values <- seq(0.2, 1.0, by = 0.01)
h1_values=0.1
h2_values=0.1

GK1<-matrix(0,M,N)  
GK2<-matrix(0,M,N) 
GK3<-matrix(0,M,N)
f1<-matrix(0,M,N) 
f2<-matrix(0,M,N) 
f3<-matrix(0,M,N)

for(pp in 1:N){
  #生成n个服从正态分布的随机数
  X <- rnorm(n, mean=1, sd=1)
  epsilon <- rnorm(n, mean=0, sd=1) #epsilon
  Y<- d*exp(X)+sin(pi*X)+c*epsilon #epsilon的系数c:0-1内取
  u <- matrix(rnorm(n, mean = 0, sd = sd0), ncol = 1) #偏差，标准差sd0
  W<- X+u
  
  #设置窗宽
  hh=besth2(x0, sd0, h1_values, h2_values, Y, W, y, c, d)
  h1=hh$best_h1
  h2=hh$best_h2
  
  print(c(pp,h1,h2))

  a1=ker((x0-W)/h1)
  a2=matrix(0,n,1)
  for (j in 1:n){
    a2[j]=Kn1((x0-W[j])/h1,sd0,h1)
  }
  G1=sum(a1)
  G2=sum(a2) 
  
  m0=sum(ker((x0-W)/h1))
  m1=sum(ker((x0-W)/h1)*(W-x0))
  m2=sum(ker((x0-W)/h1)*(W-x0)^2)
  a3=m0*m2-m1^2
  G3=sum(a3)
  
  #print(c(G1,G2))
  for(xx in 1:M){
    K=ker((Y-y[xx])/h2)/h2
    GK1[xx,pp]=t(a1)%*%K
    GK2[xx,pp]=t(a2)%*%K
    GK3[xx,pp]=t(ker((x0-W)/h1)*(m2-(W-x0)*m1))%*%K
    f1[xx,pp]=GK1[xx,pp]/G1
    f2[xx,pp]=GK2[xx,pp]/G2
    f3[xx,pp]=GK3[xx,pp]/G3
  }
}

#实际的条件概率密度函数
f00=1/(c*sqrt(2*pi))*exp(-(y-d*exp(x0)-sin(pi*x0))^2/(2*c^2)) #s[]系数为1

#若含无效NaN列，删除该列
f1 <- f1[, colSums(is.nan(f1)) == 0]
f2 <- f2[, colSums(is.nan(f2)) == 0]
f3 <- f3[, colSums(is.nan(f3)) == 0]

#计算偏差RMSEs
r1=sqrt(mean((f1-f00)^2))
r2=sqrt(mean((f2-f00)^2))
r3=sqrt(mean((f3-f00)^2))

f01=rowMeans(f1[1:M, ])
f02=rowMeans(f2[1:M, ])
f03=rowMeans(f3[1:M, ])

# 将计算结果存储为数据框
results <- data.frame(y=y, f00=f00, f01=f01, f02=f02, f03=f03)

# 绘制散点图或曲线图，并设置图例颜色和标签
p <- ggplot(results, aes(x = y)) +
  geom_line(aes(y = f00, color = "True", linetype = "True"), size = 0.7) +
  geom_line(aes(y = f02, color = "DC", linetype = "DC"), size = 0.7) +
  geom_line(aes(y = f01, color = "NW", linetype = "NW"), size = 0.7) +
  geom_line(aes(y = f03, color = "LL", linetype = "LL"), size = 0.7) +
  scale_color_manual(values = c("True" = "black", "DC" = "red", 
                                "NW" = "blue", "LL" = "green"),
                     breaks = c("True", "DC", "NW", "LL")) +
  scale_linetype_manual(values = c("True" = "solid", "DC" = "dashed", 
                                   "NW" = "dotted", "LL" = "dotdash"),
                        breaks = c("True", "DC", "NW", "LL")) +
  labs(x = "y", y = TeX("$f_{Y|X}(y|x)$"), color = " ", linetype = " ") +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1), # 添加边框
    panel.grid.major = element_blank(), # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    legend.position = c(0.1, 0.85), # 设置图例位置为右上角
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key.width = unit(0.8, "cm"),  # 设置图例框的宽度为1厘米
    legend.key.height = unit(0.5, "cm"),  # 设置图例框的高度为0.5厘米
    legend.margin = margin(-9.0, 0.8, 0.8, 1.0)  # 设置图例的边距，上、右、下、左
  ) +
  scale_x_continuous(breaks = seq(-0.5, 3, by = 0.5))

# 显示图形
print(p)
print(c(r2,r1,r3))
time2<-Sys.time()
print(time2 - time1)
