POWER_CALCULATIONd <- function(model, n , sd0, a)
{
  library(HHG)
  library(minerva) 
  library(energy)
  library(foreach)
  library(doSNOW)
  cL <- makeCluster(60, type="SOCK")
  registerDoSNOW(cL)
  
  ITER = 388;
  nr.perm = 388;
  
  BB = foreach(iter = 1:ITER, .combine = "rbind" ) %dopar% { # .packages = c("HHG", "minerva", "energy"), .combine = "cbind"  
    
    set.seed(iter)
    
    #设置噪声     
    if (model<=2)
    {
      noise=rnorm(n)
    } 
    else{
      noise=1
    }
    #noise=rnorm(n)
    
    #define models
    if (model == 1)  #linear
    {	
      xs = rnorm(n);
      ys = a*xs + noise;
    }
    
    
    if (model == 2)    #quadratic
    {
      xs = rnorm(n)
      ys = a*xs^2 + noise;
    }
    
    xs=as.matrix(xs) #不可观测
    ys=as.matrix(ys)
    
    library(GeneralizedHyperbolic)
    #双指数
    mu=0; beta=sd0/sqrt(2); alpha=beta
    us <- rskewlap(n, mu, alpha, beta,  param = c(mu, alpha, beta))
    #可观测数据
    ws<- xs+us #n*1
    
    ys_null=ys[sample(n),]
    
    X = ws; #输入X替换为可观测的ws
    Y =ys;
    Y_NULL= ys_null;
    
    source("DC.R")
    dc=DC(X,Y,n,sd0,nr.perm,a,noise)
    # 调用DC函数并存储结果到 DC_results 中
    result_dc1 <- (dc[[1]]>dc[[4]])+0  #+0是一个偏差，可以调整
    result_nw1 <- (dc[[2]]>dc[[5]])+0 
    result_ll1 <- (dc[[3]]>dc[[6]])+0 
    
    result_dc2 <- (dc[[1]]>dc[[7]])+0  #+0是一个偏差，可以调整
    result_nw2 <- (dc[[2]]>dc[[8]])+0 
    result_ll2 <- (dc[[3]]>dc[[9]])+0
    
    result_dc3 <- (dc[[1]]>dc[[10]])+0  #+0是一个偏差，可以调整
    result_nw3 <- (dc[[2]]>dc[[11]])+0 
    result_ll3 <- (dc[[3]]>dc[[12]])+0
    
    result_dc4 <- (dc[[1]]>dc[[13]])+0  #+0是一个偏差，可以调整
    result_nw4 <- (dc[[2]]>dc[[14]])+0 
    result_ll4 <- (dc[[3]]>dc[[15]])+0
    
    result_dc5 <- (dc[[1]]>dc[[16]])+0  #+0是一个偏差，可以调整
    result_nw5 <- (dc[[2]]>dc[[17]])+0 
    result_ll5 <- (dc[[3]]>dc[[18]])+0
    
    c(result_dc1,result_nw1,result_ll1,result_dc2,result_nw2,result_ll2,result_dc3,result_nw3,result_ll3,result_dc4,result_nw4,result_ll4,result_dc5,result_nw5,result_ll5)
    ##c(result_dc)
  }

  stopCluster(cL)
 
  pdc1 = mean(BB[,1])
  pnw1 = mean(BB[,2])
  pll1 = mean(BB[,3])
  
  pdc2 = mean(BB[,4])
  pnw2 = mean(BB[,5])
  pll2 = mean(BB[,6])
  
  pdc3 = mean(BB[,7])
  pnw3 = mean(BB[,8])
  pll3 = mean(BB[,9])
  
  pdc4 = mean(BB[,10])
  pnw4 = mean(BB[,11])
  pll4 = mean(BB[,12])
  
  pdc5 = mean(BB[,13])
  pnw5 = mean(BB[,14])
  pll5 = mean(BB[,15])
  
  
  power = list(pdc = c(pdc1, pdc2, pdc3, pdc4, pdc5),
               pnw = c(pnw1, pnw2, pnw3, pnw4, pnw5),
               pll = c(pll1, pll2, pll3, pll4, pll5))
  
  print(power)
  
  return(list(power = power))

}

POWER_CALCULATIONd(2, 100 , sqrt(0.25), 0)

