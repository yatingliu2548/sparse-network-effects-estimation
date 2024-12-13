% Assuming model_test_nondeg and two_sample_estimator are already defined in MATLAB#
%n1=600;
%n2=600;
n1=100;
%rho_sizes = [0.1,0.3,0.5,0.7,0.9,1]; %n2-n1; or rho2;
rho_sizes = [0,100,200,400,600,800]; %n2-n1; or rho2;
rho1=0.9;
rho2=0.9;
%rho1=rho;
%rho2=rho;
alpha1=1.8;
alphasizes = [1.2,1.8]; %alpha2
deg_choose=[0,1,2]; % 0: both are degenerate, 1: one is degenerate, 2: both are non-degenerate
dist1=1;
type1="normal";
distributionfilter=[1,2,4];
dist_notation=[10,12,13]; %10: both are normal,12: one is binomial, 13: one is uniform
replication = 1000; 
eta_num_choose=[2,3,5];
% Number of replications
parm=7+16+8+8+1; 
% Preallocate matrices for test statistics
test_stat_1 = zeros(replication*length(rho_sizes)* length(deg_choose) * length(distributionfilter)*length(alphasizes)*length(eta_num_choose),parm);
r=3;
coln = 0;
seed=123;
rng(seed); % For reproducibility, equivalent to set.seed in R

parpool; % Start the parallel pool
%% Simulation
%xi21=[36,43,52,63,80,257]
%xi31=[0,0.005,0.08,0.5,2,50]
%xi51=[0,0.015,0.24,1.5,6,150]
sigma=0.01;
for num_idex =1:length(eta_num_choose)
    eta_num=eta_num_choose(num_idex);
    for c_idx = 1:length(deg_choose)
        c= deg_choose(c_idx);
        if c==0
            d1=1;
            c1=0;
            d2=1;
            c2=0;
            eta1=0;
            eta2=0;
            Delta=0;
            xi1_val=0;
            xi2_val=0;
        elseif c==1
            d1=0;
            c1=sqrt(5);
            d2=1;
            c2=0;
            eta1=5;
            eta2=0;
            Delta=5;
            xi2_val=0;
            if eta_num==2
                xi1_val=257;
            elseif eta_num==3    
                xi1_val=50;        
            elseif eta_num==5  
                xi1_val=150;        
            end
        elseif c==2
            d1=0;
            d2=0;
            c1=sqrt(5);
            c2=sqrt(5);
            eta1=5;
            eta2=5;
            Delta=0;
            if eta_num==2
                xi1_val=257;
                xi2_val=257;
            elseif eta_num==3    
                xi1_val=50; 
                xi2_val=50;       
            elseif eta_num==5  
                xi1_val=150;
                xi2_val=150;        
            end
        end
        for alpha_idx=1:length(alphasizes)
            alpha2 = alphasizes(alpha_idx);
    
            %city_size = city_sizes(city_idx);
            for dist_idx = 1:length(distributionfilter)
                dist2 = distributionfilter(dist_idx);
                if dist2==1
                    type2="normal";
                elseif dist2==2
                    type2="poisson";
                elseif dist2==3
                    type2="binomial";
                elseif dist2==4
                    type2="uniform";
                end
                dist_no=dist_notation(dist_idx);
                for city_idx = 1:length(rho_sizes)
                    n2 = rho_sizes(city_idx)+n1;
                    %rho2 = rho_sizes(city_idx);
                    coln = coln + 1;
        
                    results = zeros(replication,parm); % Temporarily store results for each replication
        
                    parfor i = 1:replication
                
                        %Generate error terms for two city sizes
                        [e1_ij,omega1_ij] = model_sim(n1, rho1, type1,eta_num,sigma,c1,d1);
                        [e2_ij,omega2_ij] = model_sim(n2, rho2, type2,eta_num,sigma,c2,d2);
                    
                        [hatrho1_J,hateta1_J,var1_source_Gamma1,var1_source_Gamma2,jn1,var1_source_Gamma1_th,bound1_1,bound1_2] = eta_estimation(e1_ij,omega1_ij,alpha1,r,i,eta_num);
                        [hatrho2_J,hateta2_J,var2_source_Gamma1,var2_source_Gamma2,jn2,var2_source_Gamma1_th,bound2_1,bound2_2] = eta_estimation(e2_ij,omega2_ij,alpha2,r,i,eta_num);
                        hatDelta_J=hateta1_J-hateta2_J;
                        var_total=var1_source_Gamma1+var1_source_Gamma2+var2_source_Gamma1+var2_source_Gamma2;
                        var_th=var1_source_Gamma1_th+var2_source_Gamma1_th+var1_source_Gamma2+var2_source_Gamma2;
                        results(i, :) = [eta_num,Delta,hatDelta_J,var_total,var_th,c,dist_no,...  %7
                                            n1,n2,alpha1,alpha2,c1,c2,d1,d2,rho1,rho2,eta1,eta2,xi1_val,xi2_val,dist1,dist2,... %16
                                            hatrho1_J,hateta1_J,var1_source_Gamma1,var1_source_Gamma2,jn1,var1_source_Gamma1_th,bound1_1,bound1_2,... %8
                                            hatrho2_J,hateta2_J,var2_source_Gamma1,var2_source_Gamma2,jn2,var2_source_Gamma1_th,bound2_1,bound2_2,... %8
                                            i];
                    end
                    test_stat_1((coln-1)*replication+1:(coln)*replication,:)=results;
                    fprintf('sample size %.2f alpha %f num %f distribution %f.\n', n2, alpha2,eta_num,dist_no)

            
                end
            end
        end
    end
end

delete(gcp('nocreate')); % Shut down the parallel pool



T = array2table(test_stat_1, 'VariableNames', {'eta_num','Delta','hatDelta_J','var_total','var_th','c','dist_no',...  %7
    'n1','n2','alpha1','alpha2','c1','c2','d1','d2','rho1','rho2','eta1','eta2','xi1_val','xi2_val','dist1','dist2',... %16
    'hatrho1_J','hateta1_J','var1_source_Gamma1','var1_source_Gamma2','jn1','var1_source_Gamma1_th','bound1_1','bound1_2',... %8
    'hatrho2_J','hateta2_J','var2_source_Gamma1','var2_source_Gamma2','jn2','var2_source_Gamma1_th','bound2_1','bound2_2',... %8
    'replication'});
name="eta_two_sample"+"_rho_"+rho1+"_lambda_1.2"+".csv";
writetable(T, name, 'Delimiter', ';','WriteVariableNames', true);

