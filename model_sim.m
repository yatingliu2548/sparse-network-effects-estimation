function [e_ij, omega_ij] = model_sim(n, rho, type,effect,sigma,c,d)
    % Variables and Priors
    e_ij_save = zeros(n, n);
    
    % Generate a_i and b_j % Generate Epsilon_ij
    if type =="normal"
        a_i_save = normrnd(1, 1, [n, 1]);
        gamma_ij_save = normrnd(1, 1, [n, n]);
        b_j_save = normrnd(1, 1, [n, 1]);
        epsilon_ij_save = normrnd(0, sigma, [n, n]);
    elseif type =="poisson"
        a_i_save = poissrnd(1, [n, 1]);
        b_j_save = poissrnd(1, [n, 1]);
        gamma_ij_save =  poissrnd(1, [n, n]);
        epsilon_ij_save = poissrnd(1, [n, n]);
    elseif type =="binomial"
        a_i_save =  binornd(1, rand(n, 1));
        b_j_save = binornd(1, rand(n, 1));
        gamma_ij_save =  binornd(1, rand(n, n));
        epsilon_ij_save = zeros(n, n);
    elseif type =="uniform"
        a_i_save = rand(n, 1);
        b_j_save = rand(n, 1);
        gamma_ij_save =  rand(n, n);
        epsilon_ij_save = rand(n, n);
    end
    gamma_ij_save = triu(gamma_ij_save, 1) + triu(gamma_ij_save, 1)';

    % Set the diagonal elements to zero
    gamma_ij_save(logical(eye(n))) = 0;
    
    epsilon_ij_save(logical(eye(size(epsilon_ij_save)))) = 0;

    omega_ij_save = binornd(1, rho, [n, n]);
    omega_ij_save(logical(eye(size(omega_ij_save)))) = 0;

    % Generate Error Terms
    if effect ==2
         e_ij_save =(1-d).*repmat(a_i_save, 1, n) + (1-d).*repmat(b_j_save', n, 1) + epsilon_ij_save+c.*gamma_ij_save; %eta2=c^2 and xi21=36(c=0,c=sqrt(0.05), xi21=43) if d=0; (nondeg), eta2=0 and xi21=0 if d=1(deg)
    elseif effect ==3
         e_ij_save = c.*repmat(a_i_save, 1, n)+epsilon_ij_save; %eta3=c^2, eta3=0 if c=0
        
    elseif effect==5
         e_ij_save = c.*(repmat(a_i_save, 1, n)-d).*(repmat(a_i_save.', n, 1)-d)+epsilon_ij_save;  %eta5=[0.5,1,5] if c=[0.5,1,5]^{1/2} d=0(nondeg); eta5=0; eta5=0, xi51=0 if c=1, d=mean(a_i_save)(deg) 
    end
   
    e_ij_save = e_ij_save .* omega_ij_save;
    e_ij_save(logical(eye(size(e_ij_save)))) = NaN;
    e_ij=e_ij_save;
    omega_ij=omega_ij_save;
end
