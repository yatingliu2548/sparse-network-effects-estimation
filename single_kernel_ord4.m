 function single_kernel_ord4_output = single_kernel_ord4(e_ijkl,a_n,b_n,c_n)
    h_stars_save = NaN(1, 6);
    xi_matr = e_ijkl;
    xi_matr(logical (eye(size(xi_matr)))) = NaN; % Set diagonal to NaN
    xi_matr_tr=xi_matr.';
    elements = xi_matr_tr(~isnan(xi_matr_tr));
        
   % h1_star
    h_stars_save(1) = mean(elements.^2);
    % h2_star
    tmp_1=e_ijkl;
    tmp_2=e_ijkl;
    tmp_1(triu(true(size(tmp_1)))) = NaN;
    tmp_2(tril(true(size(tmp_2)))) = NaN;
    tmp_1=tmp_1.';
    a_2_clt = tmp_1(~isnan(tmp_1));
    b_2_clt = tmp_2(~isnan(tmp_2));
    h_stars_save(2) = mean(a_2_clt .* b_2_clt);
    % h3_star
    xi_vec = elements; % Assuming elements is xi.vec
    a_3 = repmat(xi_vec,2,1);
    a_3=a_3(:);
    B_temp = reshape(xi_vec, 3, 4); % Adjust dimension if needed
    combine_temp_3 = B_temp(2:end, :); % Adjust logic based on R code's intention
    for i = 2:3
        combine_temp_3 = [combine_temp_3; B_temp([1:i-1, i+1:end], :)];  % MATLAB uses ';' for vertical concatenation (rbind equivalent)
    end
    b_3 = combine_temp_3(:);
    h_stars_save(3) = mean(a_3 .* b_3);
    
    %h4_star
    xi_vec_transp= xi_matr(~isnan(xi_matr));
    a_4 = repmat(xi_vec_transp, 2,1);
    a_4=a_4(:);
    D_temp = reshape(xi_vec_transp, 3, 4); 
    combine_temp_4 = D_temp(2:end, :);
     for i = 2:3
        combine_temp_4 = [combine_temp_4; D_temp([1:i-1, i+1:end], :)];  % MATLAB uses ';' for vertical concatenation (rbind equivalent)
    end
    b_4 = combine_temp_4(:);
    h_stars_save(4) = mean(a_4 .* b_4);
   % h5_star
   h_stars_save(5) = mean([a_3; b_4] .* [b_4; a_3]);
   %h6_star
   e=e_ijkl;
   h_stars_save(6) = (e(1,2)*e(3,4) + e(1,2)*e(4,3) + e(2,1)*e(3,4) + e(2,1)*e(4,3) + ...
                  e(1,3)*e(2,4) + e(1,3)*e(4,2) + e(3,1)*e(2,4) + e(3,1)*e(4,2) + ...
                  e(1,4)*e(2,3) + e(1,4)*e(3,2) + e(4,1)*e(2,3) + e(4,1)*e(3,2)) / 12;
    
    h_star = a_n * (h_stars_save(1) + h_stars_save(2)) + ...
                 c_n * h_stars_save(6) + ...
                 b_n * (h_stars_save(3) + h_stars_save(4) + 2 * h_stars_save(5));
                 
    % Preallocate the output vector with NaNs
    single_kernel_ord4_output = NaN(1,5);

    % Subtract h_star from each element in h_stars_save(2) to h_stars_save(5)
    % and assign to the output vector
    single_kernel_ord4_output(5) = h_stars_save(5) - h_star;
    single_kernel_ord4_output(4) = h_stars_save(4) - h_star;
    single_kernel_ord4_output(3) = h_stars_save(3) - h_star;
    single_kernel_ord4_output(2) = h_stars_save(2) - h_star;
end
    
  