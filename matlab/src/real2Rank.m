function [rnk] = real2Rank(X)
    [n, d] = size(X);
    rnk = zeros(n, d);
    for i = 1:n
        x = X(i, :);
        nan_index = isnan(x);
        [ ~, ~, r] = unique(x);
        %change it to score[0 1]
        r(nan_index) = NaN;%keep the missing point
        rnk(i, :) = (r-1)/(max(r)-1);
    end
end