# DENFIS
my_denfis<-function(trn,tst,cost_range,gamma_range,epsilon_range)
{  
  MSE.mat<-array(0, dim=c(length(cost_range), length(gamma_range), length(epsilon_range)))
  count1=1
  # with 5-fold CV
  for (cost in cost_range)
  {
    count2=1  	
    for (gamma in gamma_range)
    {
      count3=1			
      for (epsilon in epsilon_range)
      {
        m<-svm(trn$data, trn$labels, cost=cost, epsilon=epsilon, gamma=gamma, cross=5)
        # average MSE
        MSE.mat[count1,count2,count3]<-median(m$MSE)
        count3=count3+1
      }
      count2=count2+1
    }
    count1=count1+1
  }
  # select best performed m
  bestparam<-which(MSE.mat==min(MSE.mat),arr.ind=T)
  best.param<-c(cost_range[bestparam[1,1]], gamma_range[bestparam[1,2]], epsilon_range[bestparam[1,3]])
  # build best SVR
  m<-svm(trn$data, trn$labels, cost=best.param[1], gamma=best.param[2], epsilon=best.param[3])
  # test
  pred<-predict(m, tst$data)
  # error measures
  MASE<-my_mase(tst$data, tst$labels, pred)
  res<-tst$labels-pred #mean
  RMSE<-sqrt(mean(res^2))
  MAPE<-mean(abs(res/mean(tst$labels)))*100
  return(list(model=m, MSE=MSE.mat, best.param=best.param, pred=pred, res=res, MASE=MASE, RMSE=RMSE, MAPE=MAPE))
}