function h = calculate_g11_J(A, i, m,combs)
    % Validate input
    n=size(A,1);
    if i > m
        error('Index i must be less than or equal to m');
    end

    % Filter combinations to include those where i is in the combinations
    valid_combs = combs(any(combs == i, 2), :);
       % Calculate the denominator
    denominator = size(valid_combs, 1);

    % Calculate the sum of products for the filtered combinations
    sum_products= zeros(denominator,1);
    r=size(valid_combs, 2);
    doubles=nchoosek(1:r, 2);
    for t = 1:size(doubles, 1)
    
        idx_1=doubles(t,1);
        idx_2=doubles(t,2);
         
        i_2=valid_combs(:, idx_1);
        j = valid_combs(:, idx_2);
        idx_ij = i_2+(j-1)*n;
		idx_ji = j+(i_2-1)*n;
        e_ij=A(idx_ij);
        e_ji=A(idx_ji);
        sum_products = sum_products + (e_ij +e_ji)/(2*nchoosek(r,2));
  
    end

    % Compute h(X_i)
    if denominator > 0
        h = sum(sum_products) / denominator;
    else
        h = 0; % To avoid division by zero
    end
end
