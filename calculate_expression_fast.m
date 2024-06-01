function result = calculate_expression_fast(A, B)
    n = size(A, 1); % Assuming A and B are n x n matrices
    if size(B, 1) ~= n || size(A, 1) ~= size(B, 1)
        error('Matrices A and B must be of the same dimension.');
    end
    
    sum_value = 0;
    for i = 1:n
        for j = 1:n
            if i ~= j
                % Sum A_ik where k is not equal to i or j
                k_indices = setdiff(1:n, [i, j]); % Indices excluding i and j
                sum_Aik = sum(A(i, k_indices));
                
                % Calculate term for i, j and accumulate
                term = (sum_Aik / (n - 2)) * B(i, j);
                sum_value = sum_value + term;
            end
        end
    end
    
    % Normalize by the number of 2-element combinations from n elements
    denominator = nchoosek(n, 2);
    result = sum_value / denominator;
end
