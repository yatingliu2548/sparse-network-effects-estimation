% Load input data
load('input_data.mat');  % This should include e, omega, alpha, r, seed, eta_num
e=double(e);
% Call the function
[hatrho_J, hateta_J, var_source_Gamma1, var_source_Gamma2, jn, var_source_Gamma1_th, bound1, bound2, g31iJ,g31iJ_2] = eta_estimation2(e, omega, alpha, r, seed, eta_num);

% Save the outputs
save('output_data.mat', 'hatrho_J', 'hateta_J', 'var_source_Gamma1', 'var_source_Gamma2', ...
     'jn', 'var_source_Gamma1_th', 'bound1', 'bound2',"g31iJ","g31iJ_2");