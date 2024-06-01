function [numerator_5,sigma_square_51] = eta5_nondegen_concentration_test(n, error_mari)
    % Data Manipulation
    xi_matr = error_mari;
    xi_matr(logical(eye(size(xi_matr)))) = NaN; % Set diagonal to NaN
    xi_matr_tr=xi_matr';
    xi_vec = xi_matr_tr(~isnan(xi_matr_tr)); % Remove NaN values

    
    % '@Var(e_ij)_mmt_estimation
    zero_diag = error_mari;
    zero_diag(logical(eye(size(zero_diag)))) = 0; % Set diagonal to zero
    
    hat_g11 = NaN(n,1);
    for i = 1:n
        hat_g11(i) = sum(zero_diag(i,:) + zero_diag(:,i)')/(n-1)/2;
    end
    hat_g11 = hat_g11 - mean(xi_vec);
    
    % '@Cov(e_ij,e_il)=eta_3_with_count=(n^2-n)(n-2)
    a_3 = repmat(xi_vec',n-2,1);
    a_3=a_3(:);
    
    % '@Cov(e_ij,e_kj)=eta_4_with_count=(n^2-n)(n-2)
    xi_vec_transp = xi_matr(~isnan(xi_matr)); 
    D_temp = reshape(xi_vec_transp, n-1, n);
    combine_temp_4 = D_temp(2:end,:);
    for i = 2:(n-1)
        combine_temp_4 = [combine_temp_4;D_temp([1:i-1, i+1:end],:)]; % rbind in MATLAB
    end
    b_4 = combine_temp_4(:); % Flatten to vector
    
    % '@Cov(e_ij,e_ki)((e_ij,e_jk))
    a_3b_4=[a_3;b_4];
    b_4a_3=[b_4; a_3];
    numerator_5 = mean(b_4.*a_3) - mean(xi_vec)^2;
    
    hat_g51_1 = NaN(n, 1);
    hat_g51_2 = NaN(n, 1);
    hat_g51_3= NaN(n, 1);
    for i = 1:n
        xi_i=xi_matr(i, :);
        xi_i=xi_i(~isnan(xi_i));
        xi_i_or=xi_i;
        xi_i=repmat(xi_i,n-2,1);
        xi_i=xi_i(:);
        indices = setdiff(1:size(xi_matr, 1), i);
        xi_ma=xi_matr(indices,indices);
        xi_ma_tr=xi_ma';
        xi_vec_i = xi_ma_tr(~isnan(xi_ma_tr)); 
        hat_g51_2(i) = sum(xi_i .* xi_vec_i);
      
        xi_vec_i2 = xi_ma(~isnan(xi_ma)); 
        xi_i2=xi_matr(:,i);
        xi_i2=xi_i2(~isnan(xi_i2));
        xi_i2=repmat(xi_i2',n-2,1);
        xi_i2=xi_i2(:);
        hat_g51_3(i) = sum(xi_i2 .* xi_vec_i2);

        xi_i1=repmat(xi_i_or,n-1,1);
        xi_i1=xi_i1(:);
        B=reshape(xi_i1,n-1,n-1)';
        B(logical(eye(size(B)))) = NaN; % Set diagonal to NaN
        b=B(:);
        b= b(~isnan(b));
        hat_g51_1(i) = sum(xi_i2.* b);
    end


    hat_g51 = (hat_g51_1 + hat_g51_2 + hat_g51_3) / (n-1) / (n-2) / 3;
    hat_g51 = hat_g51 -mean(a_3.*b_4);
    sigma_square_51 = mean((3 * hat_g51 - 4 * mean(xi_vec) .* hat_g11).^2)/n;
    
  
end