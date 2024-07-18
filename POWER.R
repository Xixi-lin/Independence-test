rm(list=ls())
#source("powerN.R")
source("powerD.R")
source("DC.R")
source("select_h.R")

model=1
n=100
#coef=25 #25时是独立
coef=seq(1, 25, by = 1)
M=length(coef)
sd0=sqrt(0.1)
a=0
result=matrix(0,M,6)

for(i in 1:M){
 #result[i, ] <- matrix(POWER_CALCULATIONn(model, n, coef[i], sd0, a)$power,1,6) #model, n,coef
  result[i, ] <- matrix(POWER_CALCULATIONd(model, n, coef[i], sd0, a)$power,1,6) #model, n,coef
}
f1=result[ ,1]
f2=result[ ,2]
f3=result[ ,3]
f4=result[ ,4]
f5=result[ ,5]
f6=result[ ,6]
#result <- POWER_CALCULATIONn(model, n, coef, sd0, a)
#print(result)

# 将计算结果存储为数据框
results <- data.frame(coef=coef, f1=f1, f2=f2, f3=f3, f4=f4, f5=f5, f6=f6)

# 绘制散点图或曲线图，并设置图例颜色和标签
p <- ggplot(results, aes(x = coef)) +
  geom_line(aes(y = f1), linetype = "solid", size = 0.7, color = "red") +
  geom_line(aes(y = f2), linetype = "dotted", size = 0.7, color = "blue") +
  geom_line(aes(y = f3), linetype = "dashed", size = 0.7, color = "black") +
  geom_line(aes(y = f4), linetype = "dotdash", size = 0.7, color = "green") +
  geom_line(aes(y = f5), linetype = "dotdash", size = 0.7, color = "orange") +
  geom_line(aes(y = f6), linetype = "dotdash", size = 0.7, color = "purple") +
  labs(x = "noise amplitude", y = "test power") +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1), # 添加边框
    panel.grid.major = element_blank(), # 去掉主要网格线
    panel.grid.minor = element_blank()  # 去掉次要网格线
  )+
  scale_x_continuous(breaks = seq(1, 25, by = 0.01))

# 显示图形
print(p)