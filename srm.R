library(TripleR)
library(srm)
library(doParallel)

TripleR_test <- function(city_size,error_mari){
  
  net <- matrix(NA,city_size^2-city_size,3)
  
  net[,1] <- rep(c(1:city_size), each = city_size - 1)
  indx <- matrix(rep(c(1:city_size), each = city_size),city_size,city_size)
  diag(indx) <- NA
  net[,2] <- na.omit(as.vector(t(indx)))
  net[,3] <- na.omit(as.vector(t(error_mari)))
  
  dataframe_net <- as.data.frame(net)
  
  fit_TripleR <- RR(V3 ~ V1*V2, data = dataframe_net)
  return(c(fit_TripleR$varComp$estimate, fit_TripleR$varComp$se))
}


srm_test <- function(city_size,error_mari){
  
  net <- matrix(NA,city_size^2-city_size,3)
  
  net[,1] <- rep(c(1:city_size), each = city_size - 1)
  indx <- matrix(rep(c(1:city_size), each = city_size),city_size,city_size)
  diag(indx) <- NA
  net[,2] <- na.omit(as.vector(t(indx)))
  net[,3] <- na.omit(as.vector(t(error_mari)))
  
  dataframe_net <- as.data.frame(net)
  colnames(dataframe_net) <- c('Actor','Partner','values')
  
  mf <- '
  %Person
  F1@A =~ 1*values@A
  F1@P =~ 1*values@P
  values@A ~~ 0*values@A + 0*values@P
  values@P ~~ 0*values@P 
  
  %Dyad
  F1@AP =~ 1*values@AP
  F1@PA =~ 1*values@PA
  values@AP ~~ 0*values@AP + 0*values@PA
  values@PA ~~ 0*values@PA 
  '
  
  mod1 <- srm::srm(mf, data = dataframe_net, conv_par=1e-4, maxiter=20)
  return(2*(1-pnorm(abs(mod1$coef/mod1$se))))
}


st_a_alt_normal <- function(n,C,rho){
  a_i.save <- rnorm(n, 1, 1)
  
  ###
  ### Generate_Epsilon_ij
  ###
  epsilon_ij.save <- matrix(rnorm(n^2, 0, 1), n, n)  
  
  
  ###
  ### Generate Error Terms
  ###
  e_ij.save <- C * matrix(rep(a_i.save, n), n, n) +  epsilon_ij.save
  omega_ij_save <- matrix(rbinom(n^2, 1, rho), nrow = n, ncol = n)
  e_ij.save=e_ij.save*omega_ij_save
  diag(e_ij.save) <- NA
 
  ###
  ### Write Output 
  ###  
  list(e_ij.save = e_ij.save)
}  


num_cores <- 20
registerDoParallel(cores=num_cores)  


#'@_Simulation
city_size <- 100
C_list <- sqrt(c(0,0.5,5,10,20)) 
replication <- 1000
test_stat3 <- matrix(NA,replication,5*6)
test_stat5 <- matrix(NA,replication,5*6)
test_stat5_2 <- matrix(NA,replication,5*6)
n2=c(100,200,300,500,700,900)
#n2=c(0.1,0.3,0.5,0.7,0.9,1)
coln <- 0
set.seed(68)
for (city_size2 in n2){
for (c_val in C_list){
  coln <- coln + 1
  
  r3 <- foreach(i=1:replication, .combine=rbind,.export = c("st_a_alt_normal", "TripleR_test", "RR")) %dopar% {
    
    #'@_e_ij
    Error_term.example <- st_a_alt_normal(n = city_size, C = 0,rho=1)
    Error_term.example2 <- st_a_alt_normal(n = city_size2, C = c_val,rho=1)
    
    
    #'@_SRM_ANOVA_method
    triplefit1 <- TripleR_test(city_size, Error_term.example$e_ij.save)
    triplefit2 <- TripleR_test(city_size2, Error_term.example2$e_ij.save)
    est3_1=triplefit1[1]
    est3_2=triplefit2[1]
    se3_1=triplefit1[7]
    se3_2=triplefit2[7]
    if (all(!is.nan(c(est3_1, est3_2, se3_1, se3_2))) && all(!is.na(c(est3_1, est3_2, se3_1, se3_2)))) {
      
      # Calculate delta using pooled variance of se3_1 and se3_2
      delta3 <- (est3_1 - est3_2) / sqrt(se3_1^2 + se3_2^2)
      
      # Calculate p-value based on the t-distribution with city_size degrees of freedom
      p_value3 <- 2 * (1 - pt(abs(delta3), df = city_size))
      
    }else{
      p_value3=666
    }
    
    
    est5_1=triplefit1[6]
    est5_2=triplefit2[6]
    se5_1=triplefit1[12]
    se5_2=triplefit2[12]
    if (all(!is.nan(c(est5_1, est5_2, se5_1, se5_2))) && all(!is.na(c(est5_1, est5_2, se5_1, se5_2)))) {
      
      # Calculate delta using pooled variance of se3_1 and se3_2
      delta5 <- (est5_1 - est5_2) / sqrt(se5_1^2 + se5_2^2)
      
      # Calculate p-value based on the t-distribution with city_size degrees of freedom
      p_value5 <- 2 * (1 - pt(abs(delta5), df = city_size))
      
    }else{
      p_value5=666
    }
    est7_1=triplefit1[2]
    est7_2=triplefit2[2]
    se7_1=triplefit1[8]
    se7_2=triplefit2[8]
    if (all(!is.nan(c(est7_1, est7_2, se7_1, se7_2))) && all(!is.na(c(est7_1, est7_2, se7_1, se7_2)))) {
      
      # Calculate delta using pooled variance of se3_1 and se3_2
      delta7 <- (est7_1 - est7_2) / sqrt(se7_1^2 + se7_2^2)
      
      # Calculate p-value based on the t-distribution with city_size degrees of freedom
      p_value7 <- 2 * (1 - pt(abs(delta7), df = city_size))
      
    }else{
      p_value7=666
    }
    c(p_value3,p_value5,p_value7)
  }
  print(city_size2)
  test_stat3[,coln] <- r3[,1]
  test_stat5[,coln] <- r3[,2]
  test_stat5_2[,coln] <- r3[,3]
}
}
write.csv(as.data.frame(test_stat3), "test_stat3_rho.csv", row.names = FALSE)  # Include row names
write.csv(as.data.frame(test_stat5), "test_stat5_rho.csv", row.names = FALSE) 
write.csv(as.data.frame(test_stat5_2), "test_stat5_2_rho.csv", row.names = FALSE) 
          