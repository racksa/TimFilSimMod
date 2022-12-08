function [p] = estimate_period(X)
% estimate_period(X) estimates the period (in indices) of the time series
% stored in X using the autocorrelation of the data. The output is stored
% in the integer p.

% Check that X is a (row or column) vector.
is_a_vector = (1 == size(X,1)) || (1 == size(X,2));
assert(is_a_vector, 'estimate_period: Input must be a vector.');

acf = autocorr(X);

% Estimate the period of the data by locating the first local maxima of the
% autocorrelation function. If X is multi-periodic, this function will only
% estimate the shortest non-trivial period.
for n = 2:length(X)-1
    above_left_data = acf(n-1) < acf(n);
    above_right_data = acf(n+1) < acf(n);
    if (above_left_data && above_right_data)
        p = n-1; % Shift by 1 because, for example, acf(1) represents the autocorrelation at zero lag.
        break;
    end
end

end

