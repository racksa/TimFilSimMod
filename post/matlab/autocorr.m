function [acf] = autocorr(X)
% autocorr(X) calculates the autocorrelation of the time series stored in
% X. The output is stored in the vector acf, which has the same orientation as X.

% Check that X is a (row or column) vector.
is_a_vector = (1 == size(X,1)) || (1 == size(X,2));
assert(is_a_vector, 'autocorr: Input must be a vector.');

% Calculate the autocorrelation.
acf = NaN(size(X));
X = X - mean(X);
for n = 1:length(X)
    acf(n) = sum( X(1 : end+1-n) .* X(n : end));    
end
acf = acf/acf(1);

end

