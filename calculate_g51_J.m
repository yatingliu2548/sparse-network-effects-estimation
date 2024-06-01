function h = calculate_g51_J(A, i, m,combs)
    % Validate input
    n=size(A,1);
    if i > m
        error('Index i must be less than or equal to m');
    end

    % Filter combinations to include those where i is in the combinations
    valid_combs = combs(any(combs == i, 2), :);
    denominator = size(valid_combs, 1);


    % Calculate the sum of products for the filtered combinations
    sum_products = zeros(denominator,1);
   
    r=size(valid_combs, 2);
    triplets=nchoosek(1:r, 3);
    for t = 1:size(triplets, 1)
        idx_1=triplets(t,1);
        idx_2=triplets(t,2);
        idx_3=triplets(t,3);
        i_2=valid_combs(:, idx_1);
        j = valid_combs(:, idx_2);
        k = valid_combs(:, idx_3);
        idx_ij = i_2+(j-1).*n;
		idx_ji = j+(i_2-1).*n;
        idx_ik = i_2+(k-1).*n;
		idx_ki = k+(i_2-1).*n;
        idx_kj = k+(j-1).*n;
        idx_jk = j+(k-1).*n;
        
        sum_products = sum_products + (A(idx_ij) .* A(idx_jk)+A(idx_kj) .* A(idx_ji)+ A(idx_ji).* A(idx_ik)+A(idx_ki) .* A(idx_ij)+ A(idx_ik).* A(idx_kj)+A(idx_jk) .* A(idx_ki))/(6 ...
            *nchoosek(r,3));
        
    end

    % Calculate the denominator
    
    % Compute h(X_i)
    if denominator > 0
        h = sum(sum_products) / denominator;
    else
        h = 0; % To avoid division by zero
    end
end
