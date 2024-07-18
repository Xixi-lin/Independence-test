DC<-function(W,Y0,n,sd0,N,a,b){ #W=X+u N=1000 #重排次数/组
  source("select_h.R")
  
  h1_values <- seq(0.1, 1, by = 0.01) 
  h2_values <- seq(0.5, 3, by = 0.02) 
  
  #设置窗宽
  HH<-besth2(sd0, h1_values, h2_values, W, Y0, a, b) #1-model1, 2-model2
  h1=HH$best_h1 
  h2=HH$best_h2 
  #HH<-optimalh1(sd0, h1_values, h2_values, Y0, W) #model3,4: 1-N, 2-D
  #h1=HH$optimal_h1 
  #h2=HH$optimal_h2
  #h1=0.08
  #h2=1.77 #1.29
  h3=h2
  
  #DC
  f0<-rep(n)
  f1<-rep(n)
  F0<-matrix(0,n,N) 
  F1<-matrix(0,n,N) 
  F01<-rep(N)
  #NW
  f2<-rep(n)
  F2<-matrix(0,n,N) 
  F02<-rep(N)
  #LL
  f3<-rep(n)
  F3<-matrix(0,n,N) 
  F03<-rep(N)
 
  for(i in 1:n){ 
    f0[i]=fdc(Y0[i], Y0[-i], h3)
    f1[i]=fDC(sd0, W[i], Y0[i], W[-i], Y0[-i], h1, h2)/f0[i]
    f2[i]=fNW(W[i], Y0[i], W[-i], Y0[-i], h1, h2)/f0[i]
    f3[i]=fLL(W[i], Y0[i], W[-i], Y0[-i], h1, h2)/f0[i]
  }
  # 定义一个函数来处理列
  remove_invalid_values <- function(col) {
    valid_values <- col[is.finite(col) & col > 0]
    return(valid_values)
  }
  valid_pp1 <- remove_invalid_values(f1)
  f01=1/n*sum(log(valid_pp1))
  
  valid_pp2 <- remove_invalid_values(f2)
  f02=1/n*sum(log(valid_pp2))
  
  valid_pp3 <- remove_invalid_values(f3)
  f03=1/n*sum(log(valid_pp3))
  
    for(pp in 1:N){
      Y <- sample(Y0) # 对 Y 进行重新排列
      Y <- matrix(Y, ncol = 1)
      for(xx in 1:n){ 
        F0[xx,pp]=fdc(Y[xx], Y[-xx], h3)
        F1[xx,pp]=fDC(sd0, W[xx], Y[xx], W[-xx], Y[-xx], h1, h2)/F0[xx,pp]
        F2[xx,pp]=fNW(W[xx], Y[xx], W[-xx], Y[-xx], h1, h2)/F0[xx,pp]
        F3[xx,pp]=fLL(W[xx], Y[xx], W[-xx], Y[-xx], h1, h2)/F0[xx,pp]
      }
      # 对fi[,pp]列应用函数
      valid_pp01 <- remove_invalid_values(F1[, pp])
      F01[pp]=1/n*sum(log(valid_pp01))
      valid_pp02 <- remove_invalid_values(F2[, pp])
      F02[pp]=1/n*sum(log(valid_pp02))
      valid_pp03 <- remove_invalid_values(F3[, pp])
      F03[pp]=1/n*sum(log(valid_pp03))
    }
    # 计算并存储分位数
    quantiles1 <- quantile(F01, 0.99)
    quantiles2 <- quantile(F02, 0.99)
    quantiles3 <- quantile(F03, 0.99)
    
    quantiles4 <- quantile(F01, 0.97)
    quantiles5 <- quantile(F02, 0.97)
    quantiles6 <- quantile(F03, 0.97)
    
    quantiles7 <- quantile(F01, 0.95)
    quantiles8 <- quantile(F02, 0.95)
    quantiles9 <- quantile(F03, 0.95)
    
    quantiles10 <- quantile(F01, 0.93)
    quantiles11 <- quantile(F02, 0.93)
    quantiles12 <- quantile(F03, 0.93)
    
    quantiles13 <- quantile(F01, 0.90)
    quantiles14 <- quantile(F02, 0.90)
    quantiles15 <- quantile(F03, 0.90)
   
    return(list(f01, f02, f03, quantiles1, quantiles2, quantiles3, quantiles4, quantiles5, quantiles6, quantiles7, quantiles8, quantiles9, quantiles10, quantiles11, quantiles12, quantiles13, quantiles14, quantiles15))
    
}

