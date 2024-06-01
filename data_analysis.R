

library(tidyverse)
library(dplyr)
library(tidyr)
theme_set(theme_bw(base_size = 12))
library(qqplotr)
library(ggplot2)
library(readr)
eta2_param <- read_delim("C:/Users/arvinyfw/Desktop/yating/sparse-network-effects-estimation/eta2_param.csv", 
                        delim = ";", escape_double = FALSE, trim_ws = TRUE)
eta3_param <- read_delim("C:/Users/arvinyfw/Desktop/yating/sparse-network-effects-estimation/eta5_param_r3_normal_nondeg_alter_0.5.csv", 
                         delim = ";", escape_double = FALSE, trim_ws = TRUE)


#eta5_param_r4_normal_deg_null.csv
#eta3_param_r4_normal_degen_null_second.csv
#eta3_param_r4_normal_alter_nondeg_0.04.csv
#eta3_param_r4_normal_degen_null.csv
#eta3_param_r3_normal_degen_null.csv
mean=0.5
eta3_param1=eta3_param%>%
  mutate(t_wq=if_else(sqrt(var_wq)>1*sqrt(log(n)/n),(hateta3_wq-mean)/sqrt(var_wq),
                      (hateta3_wq_J-mean)/sqrt(var_wq_J)),
         t_wq_th=if_else(sqrt(var_wq)>1*sqrt(log(n)/n),1,
                      0),
         t_our=if_else(sqrt(var_source_Gamma1_full)>sqrt(log(n)/n),(hateta3_full-mean)/sqrt(var_source_Gamma1_full),
                       (hateta3_J-mean)/sqrt(var_source_Gamma1)))


eta3_param1=eta3_param%>%
  mutate(t_our=if_else(sqrt(var_source_Gamma1)>sqrt(log(n)/n),(hateta3_J-mean)/sqrt(var_source_Gamma1),
               (hateta3_J-mean)/sqrt(var_source_Gamma1_2)))
ggplot(eta3_param1 %>% 
         filter( rho %in% c("0.2","0.4","1")
                 , n %in% c("100","200","400","600"),
                 alpha %in% c(1.2,1.75,2.25)
         ),
       aes(sample=t_our))+
  stat_qq(color="dodgerblue",size=0.6) +  # Color points by group#aes(color = rho)
  stat_qq_line(color="red",linewidth=0.5)+ # Add QQ line for reference
  facet_grid(alpha~ n,scale="free") +  # Separate plots for each group
  
  labs(title = "QQ Plot by Group",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles")


eta3_param=eta3_param%>%
  group_by(n,rho,alpha)%>%
  mutate(var1=var(hateta3_J),
         var1_full=var(hateta3_full),
         mean=mean(hateta3_full),
         mean_J=mean(hateta3_J))%>%ungroup()

eta3_param=eta3_param%>%
  group_by(n,rho,alpha)%>%
  mutate(var1_wq=var(hateta3_wq_J),
         var1_J=var(hateta3_J),
         mean_eta3=mean(hateta3_full),
         mean_eta3_J=mean(hateta3_J),
         mean_var_J=mean(var_source_Gamma1),
         mean_VAR_full=mean(var_source_Gamma1_full),
         mean_wq_var_J=mean(var_wq_J))%>%ungroup()

ggplot(eta3_param %>% 
         filter( rho %in% c("0.2","0.4","1")
                 , n %in% c("100","200","400"),
                 alpha %in% c(1.75)
         ),
       aes(sample=(hateta3_wq_J-mean)/sqrt(var_wq_J)))+
    stat_qq(color="dodgerblue",size=0.6) +  # Color points by group#aes(color = rho)
    stat_qq_line(color="red",linewidth=0.5)+ # Add QQ line for reference
    facet_grid(alpha~ n,scale="free") +  # Separate plots for each group
  
    labs(title = "QQ Plot by Group",
        x = "Theoretical Quantiles",
        y = "Sample Quantiles")

ggplot(eta3_param %>% 
         filter( rho %in% c("0.2","0.4","0.6","0.8","1")
                 , n %in% c("100","200","400","600"),
                 alpha %in% c(1.75,2.5)
         ),
       aes(sample=(hateta3_J-mean)/sqrt(var_source_Gamma2)))+
  stat_qq(color="dodgerblue",size=0.6) +  # Color points by group#aes(color = rho)
  stat_qq_line(color="red",linewidth=0.5)+ # Add QQ line for reference
  facet_grid(alpha~ n,scale="free") +  # Separate plots for each group
  
  labs(title = "QQ Plot by Group",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles")
mean=0
ggplot(eta3_param %>% 
         filter( rho %in% c("0.2","0.4","1")
                 , n %in% c("100","200","400","600"),
                 alpha %in% c(1.2,1.75)
         ),
       aes(sample=(hateta3_wq_J-mean)/sqrt(var_wq_J)))+
  stat_qq(color="dodgerblue",size=0.6) +  # Color points by group#aes(color = rho)
  stat_qq_line(color="red",linewidth=0.5)+ # Add QQ line for reference
  facet_grid(alpha~ n,scale="free") +  # Separate plots for each group
  
  labs(title = "QQ Plot by Group",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles")
         
 
summarise_eta_2=eta2_param %>% 
            group_by(n,rho)%>%
            summarise(mean_sigma1 = mean(sigma_square_21),
                      mean_sigma2 = mean(sigma_square_Gamma_2),
                      mean_sigma3 = mean(sigma_square_Gamma_3),
                      mean_term1 = mean(sqrt(term_1)),
                      mean_term2 = mean(sqrt(term_2)),
                      mean_gamma1 = mean(Gamma_1),
                      mean_gamma2 = mean(Gamma_2),
                      mean_gamma3 = mean(Gamma_3))
filter_square=eta2_param %>% 
  filter(sigma_square_21<bound_1)
filter=eta2_param %>% 
  filter(sigma_square_21<bound_2)


eta2_param_2=eta2_param%>%
  group_by(n,rho)%>%
  mutate(var2_1=var(t_1),
         var2_2=var(t_2),
         var2_3=var(t_3),
         mean_gamma_1=mean(Gamma_1),
         mean_gamma_2=mean(Gamma_2),
         mean_gamma_3=mean(Gamma_3),
         mean_term_1=mean(term_1),
         mean_term_2=mean(term_2),
         mean_term_3=mean(term_3),
         var=var(hat_eta2))%>%
  ungroup()

ggplot(eta2_param_2 %>% 
         filter(rho %in% c("0.2","0.4",1)
                 ,n %in% c("50","100","200")
         ),
       aes(x=mean_term_3,y=var2_3))+
  geom_point()+ # Add QQ line for reference
 geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  facet_grid(rho~ n,scale="free") +  # Separate plots for each group
  
  labs(title = "200 samples for term 3",
       x = "mean(Var 1)",
       y = "Var 2")
