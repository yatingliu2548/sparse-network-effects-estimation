
library(tidyverse)
library(dplyr)
library(tidyr)
theme_set(theme_bw(base_size = 8))
library(qqplotr)
library(ggplot2)
library(readr)
library(ggplot2)
library(dplyr)

# Define color-blind friendly colors (Okabe-Ito palette)
color_blind_palette <-c( "#08306b", "#1f78b4","#a6cee3", "#fdbf6f", "#e31a1c")
color_label = c(
  expression(italic(n) == 100),
  expression(italic(n) == 200),
  expression(italic(n) == 400),
  expression(italic(n) == 600),
  expression(italic(n) == 800)
)

color_label = c(
  expression(italic(rho) == 0.1),
  expression(italic(rho) == 0.3),
  expression(italic(rho) == 0.5),
  expression(italic(rho) == 0.7),
  expression(italic(rho) == 1.0)
)

eta_param <- read_delim("C:/Users/arvinyfw/Desktop/yating/two-sample-effects/eta_one_sample_n_600.csv", 
                         delim = ";", escape_double = FALSE, trim_ws = TRUE)
#eta_param <- read.csv("C:/Users/arvinyfw/Desktop/yating/two-sample-effects/eta_2_normal.csv") 




eta_param=eta_param%>%
  mutate(var_source_Gamma1_th=ifelse(var_source_Gamma1>n^(-3/2)*sqrt(log(n))+max(0,hatrho_J^(-2))*n^(-1/2-alpha)*(log(n))^(1/2),var_source_Gamma1,0),
    dist_deg_alpha=paste0(distribution,paste0(eta,alpha)),
    deg_alpha=paste0(eta,alpha))

# Set parameters for the sample sizes and lambda values
sample_sizes <- c(100,200,400,600,800)
lambda_values <- c(1.2,1.4,1.6,1.8)

rho_choose=0.5
n_choose=800
(rho_choose)^2*n_choose^(lambda_values-1)
(rho_choose)^2*n_choose^(lambda_values-3/2)
(rho_choose)^2*sample_sizes

#Gamma2 dominates Gamma1 when rho^2 n^(lambda-1)-->0
#Gamma2 dominates error of Gamma1 when rho^2 n^(lambda-3/2)-->0


######varying rho and fix n
data <- eta_param %>%
  filter(distribution %in% c(1)) #%>%  # Filter rows based on conditions

data <- data %>% group_by(rho, alpha, eta_num) %>%               # Group by n, alpha, and eta_num
  mutate(mean_eta = mean(hateta_J, na.rm = TRUE)) %>%  # Calculate mean of hateta_J for each group
  ungroup()

library(ggplot2)

# Define a list of labels for as_labeller using expression

custom_labels <- c(
  `01.2` = "lambda==1.2*','~xi[1]^2==0",
  `01.8` = "lambda==1.8*','~xi[1]^2==0",
  `51.2` = "lambda==1.2*','~xi[1]^2>0",
  `51.8` = "lambda==1.8*','~xi[1]^2>0",
  `2` = "eta[2]",
  `3` = "eta[3]",
  `5` = "eta[5]"
)


data$eta_num <- factor(data$eta_num, levels = c("2", "3", "5"))
data$deg_alpha <- factor(data$deg_alpha, levels = c("01.2", "01.8", "51.2", "51.8"))
# Create a labeller using as_labeller with the custom_labels
custom_labeller <- as_labeller(custom_labels, default = label_parsed)

# Create the QQ plot using ggplot2
ggplot(data, aes(sample = (hateta_J-eta) / sqrt(var_source_Gamma1_th + var_source_Gamma2), color = as.factor(rho))) +
  stat_qq(size = 0.3, alpha = 0.9) +  # Plot QQ points
  stat_qq_line(color = "black", linewidth = 0.6, linetype = "dashed",alpha = 0.8) +  # Add QQ line for reference
  facet_grid((eta_num) ~ (deg_alpha), scales = "free", labeller = as_labeller(custom_labeller)) +  # Apply custom labels
  labs(
    x = "Theoretical Quantiles",
    y = "Sample Quantiles",
    color = "Sparsity Value"  # Legend title for sample size
  ) +
  scale_color_manual(values = color_blind_palette, labels = color_label) +  # Use colorblind-friendly palette
  guides(color = guide_legend(override.aes = list(size = 1.5))) +  # Smaller dot size in legend
  theme_minimal() +
  theme(
    legend.position = "right",  # Adjust the legend position as needed
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10, face = "bold"),  
    axis.title = element_text(size = 10),  # Adjust axis title size
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 10, face = "bold"),  ## Adjust axis tick label size# Adjust legend text size for clarity
  )





######varying n and fix rho
data <- eta_param %>%
  filter(distribution %in% c(1)) #%>%  # Filter rows based on conditions
  
data <- data %>% group_by(n, alpha, eta_num) %>%               # Group by n, alpha, and eta_num
  mutate(mean_eta = mean(hateta_J, na.rm = TRUE)) %>%  # Calculate mean of hateta_J for each group
  ungroup()

library(ggplot2)

# Define a list of labels for as_labeller using expression

custom_labels <- c(
  `01.2` = "lambda==1.2*','~xi[1]^2==0",
  `01.8` = "lambda==1.8*','~xi[1]^2==0",
  `51.2` = "lambda==1.2*','~xi[1]^2>0",
  `51.8` = "lambda==1.8*','~xi[1]^2>0",
  `2` = "eta[2]",
  `3` = "eta[3]",
  `5` = "eta[5]"
)


data$eta_num <- factor(data$eta_num, levels = c("2", "3", "5"))
data$deg_alpha <- factor(data$deg_alpha, levels = c("01.2", "01.8", "51.2", "51.8"))
# Create a labeller using as_labeller with the custom_labels
custom_labeller <- as_labeller(custom_labels, default = label_parsed)

# Create the QQ plot using ggplot2
ggplot(data, aes(sample = (hateta_J-eta) / sqrt(var_source_Gamma1_th + var_source_Gamma2), color = as.factor(n))) +
  stat_qq(size = 0.3, alpha = 0.9) +  # Plot QQ points
  stat_qq_line(color = "black", linewidth = 0.6, linetype = "dashed",alpha = 0.8) +  # Add QQ line for reference
  facet_grid((eta_num) ~ (deg_alpha), scales = "free", labeller = as_labeller(custom_labeller)) +  # Apply custom labels
  labs(
    x = "Theoretical Quantiles",
    y = "Sample Quantiles",
    color = "Network Size"  # Legend title for sample size
  ) +
  scale_color_manual(values = color_blind_palette, labels = color_label) +  # Use colorblind-friendly palette
  guides(color = guide_legend(override.aes = list(size = 1.5))) +  # Smaller dot size in legend
  theme_minimal() +
  theme(
    legend.position = "right",  # Adjust the legend position as needed
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10, face = "bold"),  
    axis.title = element_text(size = 10),  # Adjust axis title size
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 10, face = "bold"),  ## Adjust axis tick label size# Adjust legend text size for clarity
  )










#############################################################################
##two sample
eta_param <- read_delim("C:/Users/arvinyfw/Desktop/yating/two-sample-effects/eta_two_sample_rho_0.5_.csv", 
                        delim = ";", escape_double = FALSE, trim_ws = TRUE)
#eta_param <- read.csv("C:/Users/arvinyfw/Desktop/yating/two-sample-effects/eta_2_normal.csv") 

rho1_val=1
eta_param=eta_param%>%
  mutate(deg_label=case_when(
    c == 0 ~ "b_deg",   # Both are degenerate
    c == 1 ~ "o_deg",  # One is degenerate
    c == 2 ~ "b_nondeg"    # Both are non-degenerate
  ), 
  diff_rho=case_when(
    rho2==0.1 ~ 0.1/rho1_val,
    rho2==0.3 ~ 0.3/rho1_val,
    rho2==0.5 ~ 0.5/rho1_val,
    rho2==0.7 ~ 0.7/rho1_val,
    rho2==0.9 ~0.9/rho1_val,
    rho2==1 ~ 1/rho1_val,
  ))

eta_param=eta_param%>%
  mutate(deg_label=case_when(
    c == 0 ~ "b_deg",   # Both are degenerate
    c == 1 ~ "o_deg",  # One is degenerate
    c == 2 ~ "b_nondeg"    # Both are non-degenerate
  ), 
  diff_n=case_when(
    n2==100 ~ 0,
    n2==200 ~ 100,
    n2==300 ~ 200,
    n2==500 ~ 400,
    n2==700 ~ 600,
    n2==900 ~ 800
  ),
  var1_source_Gamma1_th=ifelse(var1_source_Gamma1>n1^(-3/2)*sqrt(log(n1))+max(0,hatrho1_J^(-2))*n1^(-1/2-alpha1)*(log(n1))^(1/2),var1_source_Gamma1,0),
  var2_source_Gamma1_th=ifelse(var2_source_Gamma1>n2^(-3/2)*sqrt(log(n2))+max(0,hatrho2_J^(-2))*n2^(-1/2-alpha2)*(log(n2))^(1/2),var2_source_Gamma1,0),
  var_th=var1_source_Gamma1_th+var2_source_Gamma1_th+var1_source_Gamma2+var2_source_Gamma2)


legend_order <- c( "0","100","200", "400","600","800")

color_blind_palette <-c("#33a02c", "#08306b", "#1f78b4","#a6cee3", "#fdbf6f", "#e31a1c")
color_label = c(
  expression(italic(n[2])-italic(n[1])==0),
  expression(italic(n[2])-italic(n[1]) == 100),
  expression(italic(n[2])-italic(n[1])== 200),
  expression(italic(n[2])-italic(n[1]) == 400),
  expression(italic(n[2])-italic(n[1]) == 600),
  expression(italic(n[2])-italic(n[1])== 800)
)
color_blind_palette <-c("#33a02c", "#08306b", "#1f78b4","#a6cee3", "#fdbf6f", "#e31a1c")

legend_order <- c( "-0.4","-0.2","0", "0.2","0.4","0.5")
color_label = c(
  expression(italic(rho[2])-italic(rho[1])==-0.4),
  expression(italic(rho[2])-italic(rho[1]) == -0.2),
  expression(italic(rho[2])-italic(rho[1])== 0),
  expression(italic(rho[2])-italic(rho[1]) == 0.2),
  expression(italic(rho[2])-italic(rho[1]) == 0.4),
  expression(italic(rho[2])-italic(rho[1])== 0.5)
)


custom_labels_dist <- c(
  `10` = "Both are normal",
  `12` = "One is binomial",
  `13` = "One is uniform"
)


custom_labels_lambda <- c(
  `1.2` = "lambda[1]==1.8*','~lambda[2]==1.2",
  `1.8` = "lambda[1]==1.8*','~lambda[2]==1.8"
)



  
# Define custom labels with expressions as character strings for parsing
custom_labels <- c(
  `b_deg` = "Both are deg",
  `o_deg` = "One is deg",
  `b_nondeg` = "Both are nondeg"
)

custom_labels_eta <- c(  `2` = "eta[2]",
                         `3` = "eta[3]",
                         `5` = "eta[5]")

# Create a labeller with `as_labeller` for plain text labels and `label_parsed` for expressions
custom_labeller <- labeller(
  deg_label = as_labeller(custom_labels),  # For non-math expressions
  eta_num = as_labeller(custom_labels_eta,default=label_parsed)  # For math expressions
)

custom_labeller <- labeller(
  deg_label = as_labeller(custom_labels),  # For non-math expressions
  dist_no = as_labeller(custom_labels_dist)  # For math expressions
)

custom_labeller <- labeller(
  alpha2 = as_labeller(custom_labels_lambda,default = label_parsed),  # For non-math expressions
  eta_num = as_labeller(custom_labels_eta,default=label_parsed)  # For math expressions
)


data <- eta_param %>% filter(dist_no==10, deg_label=="b_deg")
data <- data %>% group_by(diff_rho,alpha2,eta_num)%>%
  mutate(mean_delta=mean(hatDelta_J))%>%ungroup()

# Create the QQ plot using ggplot2
ggplot(data, aes(sample = (hatDelta_J - mean_delta) / sqrt(var_th), color = as.factor(diff_rho))) +
  stat_qq(size = 0.3, alpha = 0.9) +  # Plot QQ points
  stat_qq_line(color = "black", linewidth = 0.6, linetype = "dashed", alpha = 0.8) +  # Add QQ line for reference
  facet_grid( eta_num ~alpha2, scales = "free", labeller = custom_labeller) +  # Apply custom labels
  labs(
    x = "Theoretical Quantiles",
    y = "Sample Quantiles",
    color = "Difference of Sparsity"  # Legend title
  ) +
  scale_color_manual(values = color_blind_palette,breaks=legend_order, labels = color_label) +  # Use colorblind-friendly palette
  guides(color = guide_legend(override.aes = list(size = 1.5))) +  # Smaller dot size in legend
  theme_minimal() +
  theme(
    legend.position = "right",  # Adjust legend position as needed
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8,face="bold"),
    axis.title = element_text(size = 10),  # Adjust axis title size
    axis.text = element_text(size = 10),
    legend.text.align = 0,  #
    strip.text = element_text(size = 8)  # Adjust facet label size
  )





















#############type I error and Type II error

set.seed(123)
eta_param <- read_delim("C:/Users/yichg/yating/sparse-network-effects-estimation/data/eta_two_sample_rho_0.9_lambda_1.2.csv", 
                        delim = ";", escape_double = FALSE, trim_ws = TRUE)

#eta_two_sample_n_600_rho_1_lambda_1.2.csv
#eta_two_sample_rho_0.9_lambda_1.2.csv
test_stat <- eta_param %>%
  filter(dist_no == 10, deg_label == "b_deg") %>%  # Filter rows based on conditions
  group_by(diff_n, alpha2, eta_num) %>%            # Group by diff_n, alpha2, and eta_num
  summarise(
    type_I = sum(abs(hatDelta_J/sqrt(var_th)) > qnorm(0.975, 0, 1)) / 1000  # Calculate Type I error rate
    #,count = n()  # Count occurrences where condition is met
  ) 

test_stat <- eta_param %>%
  filter(dist_no == 10, deg_label == "o_deg") %>%  # Filter rows based on conditions
  group_by(diff_n, alpha2, eta_num) %>%            # Group by diff_n, alpha2, and eta_num
  summarise(
    type_II = sum(abs((hatDelta_J-mean(hatDelta_J))/sqrt(var_th)) >=qnorm(0.975, 0, 1)) / 1000  # Calculate Type II error rate
    #,count = n()  # Count occurrences where condition is met
  ) 

# Define colorblind-friendly colors
color_blind_palette <- c("#88CCEE", "#CC6677", "#DDCC77")
color_label = c(
  expression(italic(eta[2* "," *A])-italic(eta[2* "," *B])==0),
  expression(italic(eta[3* "," *A])-italic(eta[3* "," *B])==0),
  expression(italic(eta[5* "," *A])-italic(eta[5* "," *B])==0)
)
legend_order=c(2,3,5)

lambda_label <- c(
  expression(italic(lambda[A])==1.2~ "," ~italic(lambda[B])==1.2),
  expression(italic(lambda[A])==1.2~ "," ~italic(lambda[B])==1.8)
)
lambda_legend=c(1.2,1.8)

# Assuming your dataframe is called `data`
ggplot(test_stat, aes(x = diff_n, y = 1-type_II, color = as.factor(eta_num), shape = as.factor(alpha2), group = interaction(eta_num, alpha2))) +
  geom_point(size = 2.5) +                # Use points with different shapes and colors
  geom_line(linewidth=0.7) +                         # Connect points with lines for each combination of `eta_num` and `alpha2`
  geom_hline(yintercept = 1, linetype = "dashed", color = "black",linewidth=0.8) +  # Add horizontal line at type_I = 0.05                       # Connect points with lines for each combination of `eta_num` and `alpha2`
  scale_shape_manual(values = c(16, 2),breaks=lambda_legend,labels = lambda_label) +  # Set specific shapes (16 = circle, 17 = triangle)
  scale_color_manual(values = color_blind_palette,breaks=legend_order, labels = color_label) +  # Use colorblind-friendly palette
  labs(
    x =  expression("Difference in Network Size" ~ (n[2] -n[1])),
    y = "Power",
    color = "Difference in Network Effects",
    shape =  "Difference in Subsampling"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",  # Adjust legend position as needed
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8,face="bold"),
    axis.title = element_text(size = 8),  # Adjust axis title size
    axis.text = element_text(size = 8),
    legend.text.align = 0,  #
    strip.text = element_text(size = 8)
  )

