function P = previous_l(matrix)
    alpha = 0.1;
    [~,n] = size(matrix);
    I = eye(n);

    P = matrix ./ (repmat(sum(matrix,2),1,n));

    P = (P + P') / 2;
end