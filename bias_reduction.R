
## ==================================================================================================== ##
## ==================================================================================================== ##
## === augmented_est is the function for Obtaining the Augmented Estimator of ========================= ##
## === for EHR-based Association Studies Accounting for Differential Misclassification ================ ##
## ==================================================================================================== ##
## ==================================================================================================== ##

## ========================================== ##
## ===== n: size of full data =============== ##
## ===== m: size of validation set ========== ##
## ===== x: predictors ====================== ##
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## !!!!!! categorical predictors have to be transferred into dummy variables !!!!!!
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

## ===== s: surrogate phenotypes ============ ##
## ===== y: true phenotypes (with NA) ======= ##
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## !!!!!! other than validation set, the value of y should be "NA" !!!!!!!!!!!!!!!!
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


augmented_est <- function(n, m, x, s, y){
  
  id <- which(y != 'NA')
  x <- cbind(1, x)
  
  nr <- dim(x)[2]
  rho <- m/n
  
  val_x = x[id,]           # validation set of x
  val_y = y[id]           # validation set of y
  val_s = s[id]           # validation set of s
  
  # glm for validation set only with true disease
  glm1 <- glm(val_y~val_x[,-1], family="binomial")
  beta_hat <- glm1$coefficients
  var_beta_hat <- summary(glm1)$coefficients[,2]**2
  
  # glm for validation set only with surrogate
  glm2 <- glm(val_s~val_x[,-1],family="binomial")
  gamma_hat <- glm2$coefficients
  var_gamma_hat <- summary(glm2)$coefficients[,2]**2
  
  # glm for all surrogate
  glm3 <- glm(s~x[,-1],family="binomial")
  gamma_bar <- glm3$coefficients
  var_gamma_bar <- summary(glm3)$coefficients[,2]**2
  
  ########################################################
  ##### Plugin estimator
  t <- table(val_s,val_y) # only validation set of y is obvered, 
  # so we can only use y in the validation set to do the estimation
  est_pooled_a1 <- t[2,2]/sum(t[,2])   # estimate sensitivity
  est_pooled_a2 <- t[1,1]/sum(t[,1])   # estimate specificity
  
  nloglik=function(t)
  {
    beta = t[1:nr];
    alpha1 = t[length(t)-1]; alpha2 = t[length(t)]
    tem= beta %*% t(x)
    p = (1-alpha2)+(alpha1+alpha2-1)*exp(tem)/(1+exp(tem))
    likelihood= s*log(p)+(1-s)*log(1-p)
    return(-sum(likelihood))
  }
  
  beta_plugin_temp <- nlminb(c(beta_hat,est_pooled_a1,est_pooled_a2),
                             nloglik,lower=c(rep(-5,nr),est_pooled_a1,est_pooled_a2),
                             upper=c(rep(5,nr),est_pooled_a1,est_pooled_a2))
  
  beta_plugin <-beta_plugin_temp$par[1:nr]
  
  Zp <- beta_plugin %*% t(x)
  
  info_matrix <- matrix(0 , nrow = nr, ncol = nr)     # information matrix
  for (i in 1:nr){
    for (j in 1:nr){
      if (i == j){
        info_matrix[i,j] <- sum( x[,i]^2 * exp(Zp) / (1+exp(Zp))^2 ) /n
      }else
        info_matrix[i,j] <- sum( x[,i] * x[,j] * exp(Zp) / (1+exp(Zp))^2) /n
      info_matrix[j,i] <- info_matrix[i,j]
    }
  }
  
  var_beta_plugin_matrix <- inv(info_matrix)/n
  var_beta_plugin <- diag(var_beta_plugin_matrix)
  ####################################################
  
  # to get the augmented estimator
  Z1 <- beta_hat %*% t(x)
  Z2 <- gamma_hat %*% t(x)
  Z3 <- gamma_bar %*% t(x)
  val_Z1 <- Z1[id]
  val_Z2 <- Z2[id]
  val_Z3 <- Z3[id]
  
  
  
  info_matrix1 <- matrix(0 , nrow = nr, ncol = nr)     # information matrix
  for (i in 1:nr){
    for (j in 1:nr){
      if (i == j){
        info_matrix1[i,j] <- sum( val_x[,i]^2 * exp(val_Z1) / (1+exp(val_Z1))^2 ) /m
      }else
        info_matrix1[i,j] <- sum( val_x[,i] * val_x[,j] * exp(val_Z1) / (1+exp(val_Z1))^2) /m
      info_matrix1[j,i] <- info_matrix1[i,j]
    }
  }
  
  info_matrix2 <- matrix(0 , nrow = nr, ncol = nr)     # information matrix
  for (i in 1:nr){
    for (j in 1:nr){
      if (i == j){
        info_matrix2[i,j] <- sum(x[,i]^2 * exp(Z3) / (1+exp(Z3))^2 ) /n
      }else
        info_matrix2[i,j] <- sum(x[,i] * x[,j] * exp(Z3) / (1+exp(Z3))^2) /n
      info_matrix2[j,i] <- info_matrix2[i,j]
    }
  }
  
  p1 <- matrix(0, nrow = nr, ncol =nr)
  p2 <- matrix(0, nrow = nr, ncol =nr)
  p3 <- matrix(0, nrow = nr, ncol =nr)
  
  for (i in 1:m){
    
    phi1 <- (val_y[i] - (exp(val_Z1[i])/(1+exp(val_Z1[i])))) %*% as.numeric(val_x[i,])
    phi1 <- inv(info_matrix1) %*% as.vector(phi1)
    p1 <- p1 + phi1 %*% t(phi1)
    
    phi2 <- (val_s[i] - (exp(val_Z3[i])/(1+exp(val_Z3[i])))) %*% as.numeric(val_x[i,])
    phi2 <- inv(info_matrix2) %*% as.vector(phi2)
    p2 <- p2 + phi2 %*% t(phi2)
    
    p3 <- p3 + phi1 %*% t(phi2)
    
  }
  
  sigma <- p1 / m
  sigma_star <- (1-rho) * p2 / m
  omega <- (1-rho) * p3 / m
  
  beta_aug = t(beta_hat - omega %*% inv(sigma_star) %*% (gamma_hat - gamma_bar))
  var_beta_aug = (sigma - omega %*% inv(sigma_star) %*% t(omega)) /m
  var_beta_aug = diag(var_beta_aug)
  
  #estimators
  est <- rbind(beta_hat, gamma_bar, beta_plugin, beta_aug)
  
  #confidence interval
  betahat.ci <- c(t(beta_hat + qnorm(0.975)*sqrt(var_beta_hat)%*%t(c(-1,1))))
  gammabar.ci <- c(t(gamma_bar + qnorm(0.975)*sqrt(var_gamma_bar)%*%t(c(-1,1))))
  beta_plugin.ci <- c(t(beta_plugin + qnorm(0.975)*sqrt(var_beta_plugin)%*%t(c(-1,1))))
  betaaug.ci <- c(t(c(beta_aug) + qnorm(0.975)*sqrt(var_beta_aug)%*%t(c(-1,1))))
  
  ci <- rbind(betahat.ci, gammabar.ci, beta_plugin.ci, betaaug.ci)
  
  # variance
  var <- rbind(var_beta_hat, var_gamma_bar, var_beta_plugin, var_beta_aug)
  
  return(list(EST = est, CI = ci, VAR=var))
}

## ===================================================================================== ##
## ===================================== Example ======================================= ##
## ===================================================================================== ##
n = 3000
m = 400

beta <- c(-0.5, 0.25, 0.38, 0.11, 1)
x1 = rbinom(n,1,0.4) # exposure group
x2_o = rnorm(n, 30, 5) # age
x3_o = rnorm(n, 150,3) # wegiht
x4 = rbinom(n,1,0.8) # race

x2 <- x2_o - mean(x2_o)
x3 <- x3_o - mean(x3_o)

x = cbind(x1,x2,x3,x4)
z = beta %*% t(cbind(1,x))
pr = 1/(1+exp(-z))     # pass through an inv-logit function
y = rbinom(n,1,pr)     # bernoulli response variable

alpha1 = 0.85          # sensitivity pr(s=1|y=1) = 0.9 for non-exposure group
alpha2 = 0.90  

a1 <- 0.90 # sensitivity 
a2 <- 0.85 # specificity 
pr_s = vector(mode = "numeric",length = n)
pr_s[x1==1] = a1*(y[x1==1]==1) + (1-a2)*(y[x1==1]==0)
pr_s[x1==0] = alpha1*(y[x1==0]==1) + (1-alpha2)*(y[x1==0]==0)
s = rbinom(n,1,pr_s)

id <- sample(n, m) 
y[-id] = NA

result <- augmented_est(n, m, x, s, y)


#### plot the result
CI <- result$CI
point_est <- result$EST
dev.off()

forest.plot <- function(i){
  ## xlim depends on the value of CI and point_est
  plot(0,0,type="n", xlab=bquote(paste(beta, .(i-1))), ylab="", yaxt="n", xaxt="n", xaxs="i", yaxs="i",ylim=c(1,7),
       xlim=c(min(CI[,(2*i-1)])-0.5,max(CI[,(2*i)])+0.5), main="")
  axis(1)
  mtext("95% Confidence Interval", outer=FALSE, line=0.4,cex=0.9)
  
  # segments for betahat
  segments(CI[1,2*i-1], c(4),CI[1,2*i], c(4), lwd=2,col="#F26A6A")
  segments(point_est[1,i], c(4)-0.1,point_est[1,i], c(4)+0.1, lwd=2,col="#F26A6A")
  segments(CI[1,2*i-1], c(4)-0.05,CI[1,2*i-1], c(4)+0.05, lwd=2,col="#F26A6A")
  segments(CI[1,2*i], c(4)-0.05,CI[1,2*i], c(4)+0.05, lwd=2,col="#F26A6A")
  
  # segments for gammabar
  segments(CI[2,2*i-1], c(3.7),CI[2,2*i], c(3.7), lwd=2, lty = 2,col="#F39759")
  segments(point_est[2,i], c(3.7)-0.1,point_est[2,i], c(3.7)+0.1, lwd=2,"#F39759")
  segments(CI[2,2*i-1], c(3.7)-0.05,CI[2,2*i-1], c(3.7)+0.05, lwd=2,"#F39759")
  segments(CI[2,2*i], c(3.7)-0.05,CI[2,2*i], c(3.7)+0.05, lwd=2,"#F39759")
  
  # segments for beta_plugin
  segments(CI[c(3),2*i-1], c(3.4),CI[c(3),2*i], c(3.4), lwd=2, lty = 3,col="#598DF3")
  segments(point_est[c(3),i], c(3.4)-0.1,point_est[c(3),i], c(3.4)+0.1, lwd=2,col="#598DF3")
  segments(CI[c(3),2*i-1], c(3.4)-0.05,CI[c(3),2*i-1], c(3.4)+0.05, lwd=2,col="#598DF3")
  segments(CI[c(3),2*i], c(3.4)-0.05,CI[c(3),2*i], c(3.4)+0.05, lwd=2,col="#598DF3")
  
  # segments for betabar
  segments(CI[c(4),2*i-1], c(3.1),CI[c(4),2*i], c(3.1), lwd=2,col="#BD74DC", lty = 4)
  segments(point_est[c(4),i], c(3.1)-0.1,point_est[c(4),i], c(3.1)+0.1, lwd=2,col="#BD74DC")
  segments(CI[c(4),2*i-1], c(3.1)-0.05,CI[c(4),2*i-1], c(3.1)+0.05, lwd=2,col="#BD74DC")
  segments(CI[c(4),2*i], c(3.1)-0.05,CI[c(4),2*i], c(3.1)+0.05, lwd=2,col="#BD74DC")
  
  # real value of the ith data (confidence interval)
  abline(v=beta[i],lty=2,col="darkgray")
  
  if (i == 1){
    legend("topleft",legend=c( expression(paste(beta[paste("V")])),
                               expression(paste(gamma[paste("F")])),
                               expression(paste(beta[paste("P")])),
                               expression(paste(beta[paste("A")])),
                               expression(paste("True value of ", beta))),
           lty=c(1,2,3,4,2),lwd=2,cex=0.8,col = c("#F26A6A","#F39759","#598DF3","#BD74DC","darkgrey"),bty="n")
  }
  
}

dim <- length(beta)
if (dim > 3){
  par(mfrow=c(round(dim/3),3))
}else{
  par(mfrow=c(1,dim))
}

for (i in c(1:dim)){
  forest.plot(i) 
}

