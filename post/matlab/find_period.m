function [p] = find_period(x)
% Calculates the period in indices of a 3D sequence.

x = x - mean(x,2);

norm_fac = sum(x.^2, 'all');

if norm_fac==0
    
    p = 1;
    
else
    
    ac = zeros(1, length(x) - 1);
    
    for n=1:(length(x) - 1)
        
        ac(n) = sum(x(:, 1:end-n) .* x(:, 1+n:end), 'all')/norm_fac;
        
    end
    
    p = 2;
    
    while (ac(p+1)>ac(p))||(ac(p-1)>ac(p))
        
        p = p + 1;
        
    end

end

% Do I want the first peak or the largest peak? When do I expect them to be
% different?

end

