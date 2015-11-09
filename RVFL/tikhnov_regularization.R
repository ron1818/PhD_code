# Tikhonov regularization
# http://en.wikipedia.org/wiki/Tikhonov_regularization

# inputs: A, x, b, a (alpha)
# outputs x_\alpha
# A = U \Sigma V^T\, \hat{x} = V D U^T b, D_{ii} = \frac{\sigma _i}{\sigma _i ^2 + \alpha ^2}

my.tikhonov.regularization<-function(A,x,b,a){
  browser()
  # 1) SVD on A
  tmp<-svd(A)
  u<-tmp$u
  s<-diag(tmp$d)
  v<-tmp$v
  
  # 1.1) check (ln(1/\lambda),||\lambda||) plot
  l=100
  s.max<-max(diag(s))
  s.min<-min(diag(s))
  lambda_l<-max(s.min,s.max*1e-16)
  lambda_1<-s.max
  lambda<-rep(NA,l)
  for (i in 1:l){
    lambda[i]<-lambda_1*(lambda_l/lambda_1)^((i-1)/(l-1))
  }
  
  # 2) calculate D_ii
  D=s/(s^2+a^2)
  
  # 3) calculate x_\alpha
  x_a<-v%*%D%*%t(u)%*%as.matrix(b)  
  
  # 4) return value
  retun(x_a)
}


