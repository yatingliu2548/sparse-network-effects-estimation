% Assuming model_test_nondeg and two_sample_estimator are already defined in MATLAB
rhosize = 0.5;
% city_sizes=[100,200,400,600,800];
city_size=600;
alphasizes = [1.1,1.3,1.5,1.7,1.9];
c_options=[0,5].^(1/2);
distributionfilter=[1,2];
replication = 1000; 
eta_num_choose=[2,3,5];
% Number of replications
parm=19; 
% Preallocate matrices for test statistics
test_stat_1 = zeros(replication*length(rhosize)* length(c_options) * length(distributionfilter)*length(alphasizes)*length(eta_num_choose),parm);
r=3;
coln = 0;
seed=68;
rng(seed); % For reproducibility, equivalent to set.seed in R

parpool; % Start the parallel pool

%% Simulation
%xi21=[36,43,52,63,80,257]
%xi31=[0,0.005,0.08,0.5,2,50]
%xi51=[0,0.015,0.24,1.5,6,150]
sigma=0.01;
for alpha_idx=1:length(alphasizes)
    alpha = alphasizes(alpha_idx);
    for num_idex =1:length(eta_num_choose)
        eta_num=eta_num_choose(num_idex);

        for c_idx = 1:length(c_options)
            c= c_options(c_idx);
            if c==0
                d=1;
                eta3=0;
                xi31_val=0;
            else
                d=0;
                eta3=5;
                if eta_num==2
                    xi31_val=257;
                elseif eta_num==3    
                    xi31_val=50;        
                elseif eta_num==5  
                    xi31_val=150;        
                end
            
            end
            %city_size = city_sizes(city_idx);
            for dist_idx = 1:length(distributionfilter)
                dist = distributionfilter(dist_idx);
                if dist==1
                    type="normal";
                elseif dist==2
                    type="poisson";
                elseif dist==3
                    type="binomial";
                elseif dist==4
                    type="uniform";
                end
                for city_idx = 1:length(rhosize)
                    rho = rhosize(city_idx);
                    coln = coln + 1;
        
                    results = zeros(replication,parm); % Temporarily store results for each replication
        
                    parfor i = 1:replication
                
                        %Generate error terms for two city sizes
                        [e_ij,omega_ij] = model_sim(city_size, rho, type,eta_num,sigma,c,d);
                    
                        [hatrho_J,hateta_J,var_source_Gamma1,var_source_Gamma2,jn,var_source_Gamma1_th,bound1,bound2] = eta_estimation(e_ij,omega_ij,alpha,r,i,eta_num);

                        results(i, :) = [eta_num,hatrho_J,hateta_J,var_source_Gamma1,var_source_Gamma2,var_source_Gamma1_th,bound1,bound2,jn,alpha,c,d,city_size,rho,sigma,eta3,xi31_val,dist,i];
                    end
                    test_stat_1((coln-1)*replication+1:(coln)*replication,:)=results;
                    fprintf('sample size %.2f alpha %f num %f distribution %f.\n', rho, alpha,eta_num,dist)

            
                end
            end
        end
    end
end

delete(gcp('nocreate')); % Shut down the parallel pool



T = array2table(test_stat_1, 'VariableNames', {'eta_num','hatrho_J','hateta_J','var_source_Gamma1','var_source_Gamma2','var_source_Gamma1_th','bound1','bound2','jn','alpha','c','deg_status','n','rho','sigma','eta','xi','distribution','replications'});
name="eta_one_sample__varyinglambda_n_600_rho_"+rhosize+".csv";
writetable(T, name, 'Delimiter', ';','WriteVariableNames', true);

