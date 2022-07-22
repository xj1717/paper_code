#############################################################
############################################################


## This is the code used for example in section 2. 

#rm(list=ls())

Sigma_Adj=function(p,N,struc,a1,a2,b,c,s) # generate covariance matrix and adjcent matrix. 
{
  # 这里推测是生成A。
  # A的取值有两种选择，一种是cov(X)，一种是$\sigma_\epsilon$，根据文章选取第二种。
  # 先求B，是$\sigma_\epsilon$
  B=matrix(0,p,p)
  for(i in 1:p)
  {
    for(j in i:p)
    {
      if (abs(i-j)<p)
      {
        #这里B是误差项的协方差矩阵
        B[i,j]=0.7^(abs(i-j)/3)
        B[j,i]=B[i,j]
      }
    }
  }
  set.seed(s)
  # 这里Pij是i类和j类之间有没有连边
  # 用均匀分布稍显复杂，使用伯努利分布生成更方便。
  # P11 = matrix(rbinom(N^2,1,a1),N,N)
  # table(runif(N^2,0,1)<a1)
  Paa=matrix(as.numeric(runif(N^2,0,1)<a1),N,N)#这里就是假设两个点同属于第一类，那么两个点之间是否有边，其实是P11
  Pbb=matrix(as.numeric(runif(N^2,0,1)<a2),N,N)#同属于第二类，其实为P22
  Pab=matrix(as.numeric(runif(N^2,0,1)<b),N,N)#属于1、2类，其实为P12=P21
  
  PAA=matrix(as.numeric(runif((2*N)^2,0,1)<a1),2*N,2*N)
  PBB=matrix(as.numeric(runif((2*N)^2,0,1)<a2),2*N,2*N)
  PAB=matrix(as.numeric(runif((2*N)^2,0,1)<b),2*N,2*N)
  ###################################################################
  # 这是怎么构造的，为啥这么构造？
  if(struc==1)
  {
    W=rbind(cbind(Paa,Pab,c*Paa,c*Pab),cbind(Pab,Pbb,c*Pab,c*Pbb),cbind(c*Paa,c*Pab,Paa,Pab),cbind(c*Pab,c*Pbb,Pab,Pbb))
  }
  if (struc==2)# 最终使用这个，结构很简单，每类100个。
  {
    W=rbind(cbind(PAA,PAB),cbind(PAB,PBB))
    
  }
  if (struc==3)
  {
    W=rbind(cbind(Paa,Pab,Paa,Pab),cbind(Pab,Pbb,Pab,Pbb),cbind(Paa,Pab,c*Paa,c*Pab),cbind(Pab,Pbb,c*Pab,c*Pbb))
    
  }
  diag(W)=0
  
  
  list(Sigma=B,Adj_mat=W)

}


library(MASS)
library(igraph)
library(e1071)

# 样本量100，向量维度5
N=100; p=5; 
# 生成邻接矩阵的次数？
num_sim=100 # given an adjcent matrix repeat X for num_sim times. 
# 类别数量？
num_block=2
# struc是临界矩阵的结构，可能对应2.3的三幅图
struc=2 #1:(2c); 2:(2a); 3:(2b)

# 这里mean_shift推测是文中u，代表mu_1向量的第一个元素的大小
mean_shift=1.6
# 这里rho推测是文中delta_0={0.05，0.1，0.3，0.5，0.7，1}
rho=0.1 #rho=1,0.05
# 这里a1就是文中a1
a1=0.5 
# 这里a2应该是文中a2
a2=0.1
# 这里c是文中的b
c=0.1 ## 0.2,0.05 
# 这里b就是b
b=seq(0,0.7,0.1)

# 随机数种子的序列
SEED=seq(from=100,to=30000,by=33)  # random seed 
# 执行算法的次数
N_simu=10   # generate random adjcent matrix for N_simu times. 

# 计算用时
timeused=matrix(nrow=N_simu,ncol=length(b))
# 聚类/分类误差
Error<-matrix(nrow=N_simu,ncol=length(b))
################# 
# 施行多次算法
for (xi in c(1:N_simu))
{
  # 每次选取一个种子
  Seed=SEED[xi]
  # 这里应该是针对（b）参数设定，对tau的不同取值施行算法。
  for(r in 1:length(b)){
    # 这里N/2，后面2*N，还是每类100个。
    S_A=Sigma_Adj(p,N/2,struc,rho*a1,rho*a2,rho*b[r],c,Seed)
    # B是误差的协方差阵，W是邻接矩阵，共200个点，40000个边（2.3中）
    B=S_A$Sigma
    W=S_A$Adj_mat
    
    #true label
    # 真实的类标是100个1，100个2（2.3中）
    if(num_block==2){
      True_label=c(rep(1,N),rep(2,N));
    }else{
      True_label=c(rep(1,N/2),rep(2,N/2),rep(3,N/2),rep(4,N/2));
    }
    
    ############### covariates
    mu1=c(mean_shift,rep(0,p-1)); mu2=-mu1;#设置μ
    mu3=c(0,mean_shift,rep(0,p-2))
    mu4=-mu3
    
    if(num_block==2){
      X1 <- mvrnorm(n=N, mu=mu1, Sigma=B)#生成N=100个向量node
      X2 <- mvrnorm(n=N, mu=mu2, Sigma=B)#mvnorm包模拟多元正态
      X=rbind(X1,X2)
    }else{
      X1 <- mvrnorm(n=N/2, mu=mu1, Sigma=B)
      X2 <- mvrnorm(n=N/2, mu=mu2, Sigma=B)
      X3 <- mvrnorm(n=N/2, mu=mu3, Sigma=B)
      X4 <- mvrnorm(n=N/2, mu=mu4, Sigma=B)
      X=rbind(X1,X2,X3,X4)
    }
    
    result=main(X,W,p,N,True_label,num_block,Iteration=FALSE)#main见另一个R代码
    Error[xi,r]=result[[1]]
    timeused[xi,r]=result[[2]]
  }
  
  
  
} 

result=list('Ave_err_NDR-KM'=apply(Error,2,mean),'Ave_timeused1'=apply(timeused,2,mean))

print(result)


