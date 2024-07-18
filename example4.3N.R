rm(list=ls())
time1<-Sys.time()
library(MASS)
library(ggplot2)
library(ggpubr)
library(latex2exp)
library(GeneralizedHyperbolic)
source("functions.R")

x0=1
sd0=sqrt(0.1) #sigma_0
n<- 200 #样本数
y=seq(-0.5, 2.5, by = 0.1)
M=length(y) #y点个数
N=100 #循环次数
h1_values <- seq(0.05, 0.6, by = 0.01)
h2_values <- seq(0.4, 1.4, by = 0.01)

GK1<-matrix(0,M,N)  
GK2<-matrix(0,M,N) 
GK3<-matrix(0,M,N) 
f1<-matrix(0,M,N) #NW
f2<-matrix(0,M,N) #DC
f3<-matrix(0,M,N) #LL

for(pp in 1:N){
  #生成n个服从正态分布的随机数
  X<- matrix(rnorm(n, mean=0, sd=1),ncol=1)
  e<- matrix(rnorm(n, mean=0, sd=0.5),ncol=1)
  Y<- X+e
  u <- matrix(rnorm(n, mean = 0, sd = sd0), ncol = 1) #偏差，标准差sd0
  #可观测数据
  W<- X+u #n*1
  
  #设置窗宽
  hh=besth3(x0, sd0, h1_values, h2_values, Y, W, y)
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
f00=1/(sqrt(2*pi)*0.5)*exp(-(y-x0)^2/(2*0.5^2))

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
  geom_line(aes(y = f00), linetype = "solid", size = 0.7, color = "black") +
  geom_line(aes(y = f01), linetype = "dotted", size = 0.7, color = "blue") +
  geom_line(aes(y = f02), linetype = "dashed", size = 0.7, color = "red") +
  geom_line(aes(y = f03), linetype = "dotdash", size = 0.7, color = "green") +
  labs(x = "y", y = TeX("$f_{Y|X}(y|x)$")) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1), # 添加边框
    panel.grid.major = element_blank(), # 去掉主要网格线
    panel.grid.minor = element_blank()  # 去掉次要网格线
  )+
  scale_x_continuous(breaks = seq(-0.5, 2.5, by = 0.5))

# 显示图形
print(p)
print(c(r2,r1,r3))
time2<-Sys.time()
print(time2 - time1)

