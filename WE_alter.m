function [e_ij, omega_ij] = WE_alter(n, rho, type,mode,c,sigma)
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
        epsilon_ij_save = poissrnd(sigma, [n, n]);
    end
    
    
    epsilon_ij_save(logical(eye(size(epsilon_ij_save)))) = 0;

    omega_ij_save = binornd(1, rho, [n, n]);
    omega_ij_save(logical(eye(size(omega_ij_save)))) = 0;

    % Generate Error Terms
    if mode ==2
         e_ij_save = repmat(a_i_save, 1, n) + repmat(b_j_save', n, 1) + epsilon_ij_save+c*gamma_ij_save;
    elseif mode ==3
         e_ij_save = epsilon_ij_save+c*repmat(a_i_save, 1, n);
        
    elseif mode==5
         e_ij_save =(repmat(a_i_save, 1, n)-1).*(repmat(a_i_save.', n, 1)-1) + epsilon_ij_save;
         % c*((repmat(a_i_save, 1, n)).*(repmat(a_i_save.', n, 1))) + epsilon_ij_save;
    end
   
    e_ij_save = e_ij_save .* omega_ij_save;
    e_ij_save(logical(eye(size(e_ij_save)))) = NaN;
    e_ij=e_ij_save;
    omega_ij=omega_ij_save;
end
