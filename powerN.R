POWER_CALCULATIONn <- function(model, n , sd0, a) 
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
  
  BB = foreach(iter = 1:ITER, .packages = c("HHG", "minerva", "energy"), .combine = "cbind") %dopar% {   
    
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
      if (model == 1)  #linear with symmetric additive noise
      {	
        xs = rnorm(n);
        ys = a*xs +noise;
      }
      
      
      if (model == 2)    #quadratic with addtive noise
      {
        xs = rnorm(n)
        ys = a*xs^2+noise;
      }
    
      
      xs=as.matrix(xs) #不可观测
      ys=as.matrix(ys)
    
      #生成n个服从正态分布的随机数
      us <- matrix(rnorm(n, mean = 0, sd = sd0), ncol = 1) #偏差，标准差sd0
      #可观测数据
      ws<- xs+us #n*1
     
      ys_null=ys[sample(n),]
      
      X = ws; #输入X替换为可观测的ws
      Y =ys;
      Y_NULL= ys_null;

      source("DC.R")
      dc=DC(X,Y,n,sd0,nr.perm,a,noise)
      # 调用DC函数并存储结果到 DC_results 中
      result_dc <- (dc[[1]]>dc[[4]])+0  #+0是一个偏差，可以调整
      result_nw <- (dc[[2]]>dc[[5]])+0 
      result_ll <- (dc[[3]]>dc[[6]])+0
      
      h_1=dc[[7]]
      h_2=dc[[8]]
      h_3=dc[[9]]
      
 #MIC  
        amic = mine(X, Y)$MIC
        amic_null=mine(X,Y_NULL)$MIC
 #HHG      
      Dx = as.matrix(dist((X)), diag = TRUE, upper = TRUE) 
      Dy = as.matrix(dist((Y)), diag = TRUE, upper = TRUE)
      hhg = hhg.test(Dx, Dy, nr.perm)
      aHHG = (hhg$perm.pval.hhg.sc <0.05)+0 
 #dcor      
      adcor=dcor(X,Y)
      adcor_null=dcor(X,Y_NULL)
      
      
      c(adcor, aHHG, amic, adcor_null, amic_null,result_dc,result_nw,result_ll,h_1,h_2,h_3)
  }
  #print(dim(BB))
  print(BB[9,])
  print(BB[10,])
  print(BB[11,])
  
  stopCluster(cL)
  bdcor=quantile(BB[4,],0.95);
  bmic=quantile(BB[5,],0.95);

  pdcor=mean(BB[1,]>bdcor)
  pHHG=mean(BB[2,])
  pmic=mean(BB[3,]>bmic)
  pdc=mean(BB[6,])
  pnw=mean(BB[7,])
  pll=mean(BB[8,])


  power=c(pdc, pnw, pll, pdcor, pHHG, pmic)
  
  # 计算经验标准差
  sd_dcor = sd(BB[1, ] > bdcor)
  sd_HHG = sd(BB[2, ])
  sd_mic = sd(BB[3, ] > bmic)
  sd_dc = sd(BB[6, ])
  sd_nw = sd(BB[7, ])
  sd_ll = sd(BB[8, ])
  
  power_sd = c(sd_dc, sd_nw, sd_ll, sd_dcor, sd_HHG, sd_mic)
  
  print(c(power, power_sd))
  
  return(list(power = power, power_sd = power_sd))
}


