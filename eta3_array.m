

taskid=str2double(getenv('SLURM_ARRAY_TASK_ID'));
jobid=str2double(getenv('SLURM_ARRAY_JOB_ID'));
seed=jobid;
disp(jobid)
rng(seed);
fprintf('seed %d', seed)
city_sizes = [200,400,600];%,200,400,600,800
rhosizes = [0.5,1];%,0.6,0.8,0.9
alphasizes = [1.75,1.8,2.0,2.5];
%city_sizes = [25,500];
%rhosizes = [0.1,1];
%replication = 100; % Number of replications
parm=27;
% Preallocate matrices for test statistics
test_stat_1 = zeros(length(city_sizes) * length(rhosizes)*length(alphasizes),parm);

coln = 0;
r=3;
%sigma=0.01;
%rng(68); % For reproducibility, equivalent to set.seed in R
c=0.2;
rhosize=0.4;
%parpool; % Start the parallel pool
mean=0.04;
for alpha_idx=1:length(alphasizes)
     alpha = alphasizes(alpha_idx);
    for city_idx = 1:length(city_sizes)
        city_size = city_sizes(city_idx);
        for rho_idx = 1:length(rhosizes)
            sigma = rhosizes(rho_idx);
            coln = coln + 1;
        
            results = zeros(1,parm); % Temporarily store results for each replication
        
            %parfor i = 1:replication
            % Generate error terms for two city sizes
                [e_ij,omega_ij] = WE_alter(city_size, rhosize,"normal",3,c,sigma);
                [e_ij2,omega_ij2] = WE_alter(city_size+200, rhosize+0.6,"normal",3,c,sigma);
                %E2 = WE_null(city_size + 100, rhosize + 0.1);
                Jn=getsubset(city_size,alpha,r);
                Jn2=getsubset(city_size+50,alpha,r);
                %Jn2=getsubset(city_size,alpha,4);
                % Compute test statistics
                [hatrho_J,hateta3_J,hateta3_full,var_source_Gamma1,var_source_Gamma2,var_source_Gamma3,jn,var_source_Gamma1_2,var_source_Gamma1_full] = eta3_estimator(e_ij,omega_ij,alpha,r,Jn);
                [hatrho_J2,hateta3_J2,hateta3_full2,var2_source_Gamma1,var2_source_Gamma2,var2_source_Gamma3,jn2,var2_source_Gamma1_2,var2_source_Gamma1_full] = eta3_estimator(e_ij2,omega_ij2,alpha,r,Jn2);

                
                hateta3_wq_J=0;
                var_wq_J=0;
                %[hateta3_wq_J,var_wq_J]=incomplete_Us_NO_debias_standz(city_size,Jn2,e_ij,3);
                bound_1 = sqrt(log(city_size)/city_size);
                bound_2 = (log(city_size)/city_size);
                fprintf('sample size %.2f alpha %f rhosize %f.\n', city_size, alpha,rhosize)
                results(1, :) = [0,hatrho_J,hateta3_J,hateta3_full,var_source_Gamma1,var_source_Gamma2,var_source_Gamma3,jn,var_source_Gamma1_2,var_source_Gamma1_full,hatrho_J2,hateta3_J2,hateta3_full2,var2_source_Gamma1,var2_source_Gamma2,var2_source_Gamma3,jn2,var2_source_Gamma1_2,var2_source_Gamma1_full,city_size,city_size+50,rhosize,rhosize+0.6,alpha,r,mean,sigma];
            %end
            test_stat_1(coln,:)=results;
        
        end
    end
end

%delete(gcp('nocreate')); % Shut down the parallel pool


name=sprintf('eta3/eta3_normal_two_sample_%d_%d.csv', jobid, taskid);
disp(name)
T = array2table(test_stat_1, 'VariableNames', {'deg_status','hatrho_J','hateta3_J','hateta3_full','var_source_Gamma1','var_source_Gamma2','var_source_Gamma3','jn','var_source_Gamma1_2','var_source_Gamma1_full','hatrho_J2','hateta3_J2','hateta3_full2','var2_source_Gamma1','var2_source_Gamma2','var2_source_Gamma3','jn2','var2_source_Gamma1_2','var2_source_Gamma1_full','n1','n2','rho1','rho2','alpha','r','mean','sigma'});
writetable(T, name, 'Delimiter', ';','WriteVariableNames', true);