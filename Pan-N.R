rm(list=ls())
time1<-Sys.time()
library(MASS)
library(ggplot2)
library(ggpubr)
library(latex2exp)
library(GeneralizedHyperbolic)
library(foreign)
library(openxlsx)
source("functions.R")

z<-read.xlsx("Pzfv.xlsx")
#set.seed(123) # 使用自己的种子值
#random_rows <- z[sample(nrow(z), 699), ]
#n=nrow(random_rows)
#Z<-as.data.frame(random_rows)[1:n,4:7]
n=nrow(z)
Z<-as.data.frame(z)[1:n,3:26] #z:3-10;f:11-18;v:19-26


x0=2
y=seq(-1, 3, by = 0.1) #-2, 2
M=length(y) #y点个数
#h1_values <- seq(0.1, 1, by = 0.1)
#h2_values <- seq(0.1, 1, by = 0.1)

GK1<-matrix(0,M,1)  
GK3<-matrix(0,M,1) 
GK4<-matrix(0,M,1) 
f1<-matrix(0,M,1) #NW
f3<-matrix(0,M,1) #LL
f4<-matrix(0,M,1) #DC_newh

Y <- Z[1:n, 17:24] #金属丰度9:16/径向速度17:24
# 创建一个函数来计算行平均值
row_average <- function(row) {
  return(mean(row, na.rm = TRUE))
}

# 使用apply函数计算每行的平均值
Y <- apply(Y, 1, row_average)
#标准化Y
Y <- (Y - mean(Y)) / sd(Y)
w <- Z[1:n, 1:8] #红移
# 使用apply函数计算每行的有效值数量
row_valid_counts <- apply(w, 1, function(row) sum(!is.na(row)))
# 计算整个数据表中的有效值总数
total_valid_count <- sum(row_valid_counts, na.rm = TRUE)
W=apply(w, 1, row_average) #行均值
# 计算每个值与其对应行均值的差的平方的和
squared_diff_sums <- rowSums((w - W)^2, na.rm = TRUE)
son=sum(squared_diff_sums) #分子
sd0=sqrt(son/(total_valid_count-n))
#标准化W
W <- (W - mean(W)) / sd(W)


#设置窗宽
#H=optimalh1(sd0, h1_values, h2_values, Y, W) #正态误差分布
#H1=H$optimal_h1
#H2=H$optimal_h2
#金属丰度
#H1=0.16
#H2=0.58
#径向速度
H1=0.18
H2=1.97

print(c(H1,H2))

a1=ker((x0-W)/H1)
a4=matrix(0,n,1)
for (j in 1:n){
  a4[j]=Kn1((x0-W[j])/H1,sd0,H1)
}
G1=sum(a1)
G4=sum(a4)

m0=sum(ker((x0-W)/H1))
m1=sum(ker((x0-W)/H1)*(W-x0))
m2=sum(ker((x0-W)/H1)*(W-x0)^2)
a3=m0*m2-m1^2
#a3=ker((x0-W)/h1)*(m2-(W-x0)*m1)
G3=sum(a3)

for(xx in 1:M){
  K=ker((Y-y[xx])/H2)/H2
  GK1[xx]=t(a1)%*%K
  GK3[xx]=t(ker((x0-W)/H1)*(m2-(W-x0)*m1))%*%K
  GK4[xx]=t(a4)%*%K
  f1[xx]=GK1[xx]/G1
  f3[xx]=GK3[xx]/G3
  f4[xx]=GK4[xx]/G4
}


#实际的条件概率密度函数
#f00=1/sqrt(2*pi*(1-rho^2)*SigY)*exp(-(rho*(x0-muX)/sqrt(SigX)-(y-muY)/sqrt(SigY))^2/(2*(1-rho^2)))

#若含无效NaN列，删除该列
f1 <- f1[, colSums(is.nan(f1)) == 0]
f3 <- f3[, colSums(is.nan(f3)) == 0]
f4 <- f4[, colSums(is.nan(f4)) == 0]


# 将计算结果存储为数据框
results <- data.frame(y=y, f1=f1, f3=f3, f4=f4)

# 绘制散点图或曲线图，并设置图例颜色和标签
p <- ggplot(results, aes(x = y)) +
  geom_line(aes(y = f4, color = "DC", linetype = "DC"), size = 0.7) +
  geom_line(aes(y = f1, color = "NW", linetype = "NW"), size = 0.7) +
  geom_line(aes(y = f3, color = "LL", linetype = "LL"), size = 0.7) +
  scale_color_manual(values = c("DC" = "red", 
                                "NW" = "blue", "LL" = "green"),
                     breaks = c("DC", "NW", "LL")) +
  scale_linetype_manual(values = c("DC" = "solid", 
                                   "NW" = "dotted", "LL" = "dotdash"),
                        breaks = c("DC", "NW", "LL")) +
  labs(x = "y", y = TeX("$f_{Y|X}(y|x)$"), color = " ", linetype = " ") +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1), # 添加边框
    panel.grid.major = element_blank(), # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    legend.position = c(0.1, 0.85), # 设置图例位置为右上角
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key.width = unit(0.8, "cm"),  # 设置图例框的宽度为0.8厘米
    legend.key.height = unit(0.5, "cm"),  # 设置图例框的高度为0.5厘米
    legend.margin = margin(-9.0, 0.8, 0.8, 1.0)  # 设置图例的边距，上、右、下、左
  )+
  scale_x_continuous(breaks = seq(-1, 3, by = 0.5))

# 显示图形
print(p)
#print(c(r2,r4,r1,r3))


time2<-Sys.time()
print(time2 - time1)
