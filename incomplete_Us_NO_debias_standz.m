function [t, varc] = incomplete_Us_NO_debias_standz(n, Jn, error_mari,mode)
    % Set diagonal of error_mari to NaN
    error_mari(logical(eye(size(error_mari)))) = NaN;

    % Define coefficients
    a_n = 1 / (n * (n-1));
    b_n = (n-2) / (n * (n-1));
    c_n = (n-2) * (n-3) / (n * (n-1));
    replctn = size(Jn,1);
    
    % Preallocate incomplete_Us
    incomplete_Us = NaN(replctn, 5);
    
    for i = 1:replctn
        % MATLAB does not have a direct equivalent to set.seed for reproducibility within a loop
        % Selection of 4 unique indices
        ijkl = Jn(i,:);
        error_mari_dim4 = error_mari(ijkl, ijkl);
        
        % Assuming single_kernel_ord4 is already defined and adapted for MATLAB
        incomplete_Us(i,:) = single_kernel_ord4(error_mari_dim4,a_n,b_n,c_n);
    end
    
    % Calculate test statistics
    means = mean(incomplete_Us, 1);
    vars = var(incomplete_Us, 0, 1)/replctn;
    t=means(mode);
    varc=vars(mode);
   
end
