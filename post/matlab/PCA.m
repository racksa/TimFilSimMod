function [A, V, percentage_weights, percentage_error_in_cumulative_soln] = PCA(X)
% Compute the PCA decomposition without centering the data -- I assume that
% if centering is relevant to the problem I will have done it beforehand.

[~, S, V] = svd(X, 0);

A = X*V;

V = V';

percentage_weights = 100 * diag(S).^2 / sum(diag(S).^2);

% The n-th entry of this vector contains the percentage reconstruction
% error if the solution is approximated using only the first n modes. The
% last entry corresponds to complete reconstruction and will always be 0.
percentage_error_in_cumulative_soln = 100 * sqrt(sum(diag(S).^2) - cumsum(diag(S).^2)) / norm(diag(S));

end

