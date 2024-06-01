% Assuming WE_null_normal and two_sample_estimator are already defined in MATLAB
city_sizes = [25,50,100,200,300,400,500];
rhosizes = [0.012,0.013,0.05,0.1,0.12,0.15,0.2,0.3, 0.4,0.5,0.7 0.8,0.9,1];
%city_sizes = [25,500];
%rhosizes = [0.1,1];
replication = 200; % Number of replications

% Preallocate matrices for test statistics
test_stat_1 = zeros(replication* length(city_sizes) * length(rhosizes),18);

coln = 0;
rng(68); % For reproducibility, equivalent to set.seed in R

parpool; % Start the parallel pool

for city_idx = 1:length(city_sizes)
    city_size = city_sizes(city_idx);
    for rho_idx = 1:length(rhosizes)
        rhosize = rhosizes(rho_idx);
        coln = coln + 1;
        
        results = zeros(replication,18); % Temporarily store results for each replication
        
        parfor i = 1:replication
            % Generate error terms for two city sizes
            [e_ij_original,e_ij,omega_ij] = WE_null(city_size, rhosize,"normal",2);
            %E2 = WE_null(city_size + 100, rhosize + 0.1);
            
            % Compute test statistics
            [t_1,t_2,t_3, hat_eta2, sigma_square_21, n,Gamma_1,Gamma_2,Gamma_3,term_1,term_2,term_3,sigma_square_Gamma_2,sigma_square_Gamma_3,hat_rho] = eta2_estimator(e_ij_original,e_ij,omega_ij,rhosize);
  
            bound_1 = sqrt(log(city_size)/city_size);
            bound_2 = (log(city_size)/city_size);
            
            results(i, :) = [t_1,t_2,t_3, hat_eta2, sigma_square_21, bound_1,bound_2,Gamma_1,Gamma_2,Gamma_3,term_1,term_2,term_3,sigma_square_Gamma_2,sigma_square_Gamma_3,city_size,rhosize,hat_rho];
        end
        test_stat_1((coln-1)*replication+1:(coln)*replication,:)=results;
        
    end
end

delete(gcp('nocreate')); % Shut down the parallel pool



T = array2table(test_stat_1, 'VariableNames', {'t_1','t_2','t_3', 'hat_eta2', 'sigma_square_21', 'bound_1','bound_2','Gamma_1','Gamma_2','Gamma_3','term_1','term_2','term_3','sigma_square_Gamma_2','sigma_square_Gamma_3','n','rho','hat_rho'});
writetable(T, "eta2_param.csv", 'Delimiter', ';','WriteVariableNames', true);

selectedRows = T.sigma_square_21 < T.bound_1;
filteredT = T(selectedRows, :);