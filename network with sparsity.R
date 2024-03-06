## data generation
WE_null_normal <- function(n,rho,seed=123){
  
  ###
  ### Variables and Priors
  ###
  e_ij.save <- matrix(0,n,n)
  
  ###
  ###'@Generate_ab_ij__Paras:a0,b0,pho_ab
  ###
  a_i.save <- rnorm(n, 1, 1)
  b_j.save <- rnorm(n, 1, 1)
  
  ###
  ###'@Generate_Epsilon_ij__Paras:a2,b2
  ###
  epsilon_ij.save <- matrix(rnorm(n^2, 0, 1), n, n)   
  diag(epsilon_ij.save) <- 0
  
  omega_ij.save <- matrix(rbinom(n^2,size=1,prob=rho),n,n)
  diag(omega_ij.save) <- 0
  
  ###
  ### Generate Error Terms
  ###
  e_ij.save <- matrix(rep(a_i.save, n), n, n) + 
    matrix(rep(b_j.save, n), n, n, byrow = TRUE) + epsilon_ij.save
  e_ij.save <- e_ij.save * omega_ij.save 
  diag(e_ij.save) <- NA
  
  
  ###
  ### Write Output 
  ###  
  return(list(e_ij.save = e_ij.save,omega_ij.save=omega_ij.save))
  
}
WE_null_Pois <- function(n,rho){
  
  ###
  ### Variables and Priors
  ###
  e_ij.save <- matrix(0,n,n)
  
  ###
  ###'@Generate_ab_ij__Paras:a0,b0,pho_ab
  ###
  a_i.save <- rpois(n, 1)
  b_j.save <- rpois(n, 1)
  
  ###
  ###'@Generate_Epsilon_ij__Paras:a2,b2
  ###
  epsilon_ij.save <- matrix(rpois(n^2, 1), n, n)    
  diag(epsilon_ij.save) <- 0
  
  
  ###
  ### Generate Error Terms
  ###
  omega_ij.save <- matrix(rbinom(n^2,size=1,prob=rho),n,n)
  diag(omega_ij.save) <- 0
  
  ###
  ### Generate Error Terms
  ###
  e_ij.save <- matrix(rep(a_i.save, n), n, n) + 
    matrix(rep(b_j.save, n), n, n, byrow = TRUE) + epsilon_ij.save
  e_ij.save <- e_ij.save * omega_ij.save 
  diag(e_ij.save) <- NA
  

  
  ###
  ### Write Output 
  ###  
  return(list(e_ij.save = e_ij.save,omega_ij.save=omega_ij.save))
  
}


hat_eta2_estimation <- function(e_ij,omega_ij){
  
  xi.matr <- e_ij
  diag(xi.matr) <- NA
  elements <- as.vector(t(xi.matr))
  xi.vec <- na.omit(elements) # hat U_1 star is its mean
  tmp_1 <- e_ij
  tmp_2 <- e_ij
  tmp_1[lower.tri(tmp_1)] <- NA # only upper triangle matrix without diagonal entries 
  tmp_2[upper.tri(tmp_2)] <- NA
  
  a_2_clt <- na.omit(as.vector(t(tmp_1)))
  b_2_clt <- na.omit(as.vector(tmp_2))
  hat_eta2_star <- mean(a_2_clt * b_2_clt) - mean(xi.vec) * mean(xi.vec) 
  
  ##hat rho
  xi.matr2 <- omega_ij
  diag(xi.matr2) <- NA
  elements2 <- as.vector(t(xi.matr2))
  xi.vec2 <- na.omit(elements2) 
  hat_rho <- mean(xi.vec2)
  hat_eta2 <- hat_rho^(-2) * hat_eta2_star
  
  return(list(hat_eta2_star=hat_eta2_star,hat_eta2=hat_eta2,hat_rho=hat_rho))
}

xi_square_21_estimation <- function(e_ij,omega_ij,hat_eta2,hat_eta2_star,hat_rho){
  xi.matr <- e_ij
  diag(xi.matr) <- NA
  elements <- as.vector(t(xi.matr))
  xi.vec <- na.omit(elements) 
  #####hat rho
  xi.matr2 <- omega_ij
  diag(xi.matr2) <- NA
  elements2 <- as.vector(t(xi.matr2))
  xi.vec2 <- na.omit(elements2) 
  
  ###'@Var(e_ij)_mmt_estimation
  ###
  zero_dig <- e_ij
  diag(zero_dig) <- 0
  zero_dig_omega <- omega_ij
  diag(zero_dig_omega) <- 0
  n=dim(e_ij)[1]
  
  hat_a_1i_star <- rep(NA,n)
  hat_a_1i_rho <- rep(NA,n)
  hat_a_1i_rho_2 <- rep(NA,n)
  for (i in 1:n) {
    hat_a_1i_star[i] <- sum(zero_dig[i,] + zero_dig[,i])/(n-1)/2    # is hat a_{1,i} star
    hat_a_1i_rho[i] <- sum(((zero_dig_omega[i,] + zero_dig_omega[,i])/2-hat_rho)^2)/(n-1)
    hat_a_1i_rho_2[i] <- sum((zero_dig_omega[i,] + zero_dig_omega[,i]))/(n-1)/2
    
    
  }
  hat_g_rho <- hat_a_1i_rho
  hat_g_rho_2 <- hat_a_1i_rho_2-hat_rho
  hat_g11 <- hat_a_1i_star - mean(xi.vec)
  
  ###
  ###'@Cov(e_ij,e_ji)=eta_2_with_count=(n^2-n)/2_(for_independent_CLT)
  ###
  xi.matr <- e_ij
  diag(xi.matr) <- NA
  tmp_1 <- xi.matr
  tmp_2 <- xi.matr
  tmp_1[lower.tri(tmp_1)] <- NA
  tmp_2[upper.tri(tmp_2)] <- NA
  tmp_1_rho <- xi.matr2
  tmp_2_rho <- xi.matr2
  tmp_1_rho[lower.tri(tmp_1_rho)] <- NA
  tmp_2_rho[upper.tri(tmp_2_rho)] <- NA
  
  a_2_clt_rho <- na.omit(as.vector(t(tmp_1_rho)))
  b_2_clt_rho <- na.omit(as.vector(tmp_2_rho))
  
  a_2_clt <- na.omit(as.vector(t(tmp_1)))
  b_2_clt <- na.omit(as.vector(tmp_2))
  hat_U_2_sqrt <- mean(a_2_clt^2 *b_2_clt ^2)
  hat_U_2 <- mean(a_2_clt *b_2_clt )
  hat_rho_2_sqrt <- sum(a_2_clt_rho*a_2_clt_rho- hat_rho^2)/(2*length(a_2_clt_rho))
  hat_a_2i_star <- rep(NA,n)
  hat_a_2i_star_rho <- rep(NA,n)
  hat_a_2i_star_rho_2 <- rep(NA,n)
  for (i in 1:n) {
    hat_a_2i_star[i] <- sum(zero_dig[i,] * zero_dig[,i])/(n-1)
    hat_a_2i_star_rho[i] <- sum((zero_dig_omega[i,] * zero_dig_omega[,i]-hat_rho^2)^2)/(n-1)
    hat_a_2i_star_rho_2[i] <- sum((zero_dig_omega[i,] * zero_dig_omega[,i]))/(n-1)
  }
  hat_g21 <- hat_a_2i_star - mean(a_2_clt * b_2_clt)
 hat_g2_rho<- hat_a_2i_star_rho
 hat_g2_rho_2<- hat_a_2i_star_rho_2-hat_rho^2
  sigma_square_21 <- hat_rho^(-4) * mean((2*hat_g21 - 4*mean(xi.vec)*hat_g11)^2) 
  
  #variance with hat rho
 
  

  sigma_square_rho_21 <- hat_rho^(-4) * mean((2*hat_g21 - 4*mean(xi.vec)*hat_g11 - 
                                                2*hat_rho^(-1) * hat_eta2_star *sqrt(hat_g_rho)
                                              -sqrt(hat_g2_rho)* hat_rho^(-2)*hat_U_2
                                  
                                              )^2)
  # is 1
  sigma_square_rho_2 <- hat_rho^(-4) * mean((2*hat_g21 - 4*mean(xi.vec)*hat_g11 - 
                                               2*hat_rho^(-1) * hat_eta2_star *  hat_g_rho_2
                                             -2*hat_g2_rho_2*hat_rho^(-2)*hat_U_2
                                        
                                              )^2)
    #4 * ((1-hat_rho)/hat_rho)* hat_eta2^2
  # is 6
  sigma_square_R_2<-  hat_rho^(-4) * mean((2*hat_g21 - 4*mean(xi.vec)*hat_g11 - 
                                                2*hat_rho^(-1) * hat_eta2_star * hat_g_rho_2
                                           -hat_g2_rho_2* hat_rho^(-2)*hat_U_2
                                           -hat_g21* mean(hat_g2_rho)
                                           )^2)
 # sigma_square_R_2 <- (hat_rho^(-2)-1) * hat_rho^(-2) * mean(hat_a_2i_star_rho)
  
  #denominator_2 <- sqrt(sigma_square_21/n)
  return(list(sigma_square_21=sigma_square_21,sigma_square_rho_21=sigma_square_rho_21,sigma_square_rho_2=sigma_square_rho_2,sigma_square_R_2=sigma_square_R_2))
  
}
## test
replication=1000
results=matrix(NA,replication, 6)
for (i in 1: replication){
 
  n_10_rho_0.1 <- WE_null_Pois(n=n, rho=rho)
  hat_eta2_result=hat_eta2_estimation(n_10_rho_0.1$e_ij.save,n_10_rho_0.1$omega_ij.save)
  hat_eta2_star = hat_eta2_result$hat_eta2_star
  hat_eta2 = hat_eta2_result$hat_eta2
  #results[i,1]=hat_eta2_result$hat_rho
  results[i,2]=hat_eta2
  results[i,3]=hat_eta2_star
  sigma_result=xi_square_21_estimation(n_10_rho_0.1$e_ij.save,n_10_rho_0.1$omega_ij.save,hat_eta2_result$hat_eta2,hat_eta2_result$hat_eta2_star,rho)
  results[i,4]=sigma_result$sigma_square_21
  results[i,5]=sigma_result$sigma_square_rho_21
  results[i,6]=sigma_result$sigma_square_R_2
  results[i,1]=sigma_result$sigma_square_rho_2
}
n=100
rho=0.03

log(n)/(n^2)
sqrt(log(n)/n)
(log(n)^(3/2))/(n^(1/2))
#data=results[,2]/(sqrt((results[,4]/n)+(results[,5]/choose(n,2))+(results[,6]/n)))
#data2=results[,2]/(sqrt((results[,4]/n)+(results[,5]/choose(n,2))))
mean=mean(na.omit(results[,2]))
data1=(results[,2])/(sqrt((results[,1]/(n))))
data2=(results[,2])/(sqrt((results[,6]/(n))))
data3=(results[,2])/(sqrt((results[,5]/(n))))
data=(results[,2])/(sqrt((results[,4]/(n))))
var=sqrt(results[,4]/n)
sum(var<sqrt(log(n)/n))

qqnorm(data)
abline(0,1, lty = 2, lwd = 1.5)
qqnorm(data2)
abline(0,1, lty = 2, lwd = 1.5)
qqnorm(data3)
abline(0,1, lty = 2, lwd = 1.5)

qqnorm(data1)
abline(0,1, lty = 2, lwd = 1.5)

