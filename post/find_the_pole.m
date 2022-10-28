function guess = find_the_pole(seg_pos, NFIL, NSEG)
% Find the choice of polar axis that minimises the error of the
% decomposition in az_polar_coeffs(...). % Takes inputs from during the
% execution of animate.m.

polar_axis = @(phi, theta) [sin(theta)*cos(phi); sin(theta)*sin(phi); cos(theta)];

fun_to_max = @(phi, theta) norm(az_polar_coeffs(polar_axis(phi, theta), seg_pos, NFIL, NSEG));

delta = 1e-5;

grad = @(phi, theta) [fun_to_max(phi+delta, theta) - fun_to_max(phi, theta); fun_to_max(phi, theta+delta) - fun_to_max(phi, theta)]/delta;

% Loop over a spiral seeding to find a good initial guess:

val = 0;

temp = zeros(3,1);

for k=1:100
    
    temp(3) = 2*(k-1)/99 - 1;
    
    r = sqrt(1 - temp(3)^2);
    
    if k==1 || k==100
        
        phi = 0;
        
    else
        
        phi = phi + 3.6/sqrt(r*NFIL);
        
    end
    
    temp(1) = r*cos(phi);
    temp(2) = r*sin(phi);
    
    new_val = norm(az_polar_coeffs(temp, seg_pos, NFIL, NSEG));
    
    if new_val > 0.9
        
        guess = temp;
        
        break;
        
    elseif new_val > val
        
        val = new_val;
        
        guess = temp;
        
    end
    
end

% Refine the guess (do it in terms of polar coordinates):

phi = atan2(guess(2), guess(1));
theta = acos(guess(3));

curr_grad = grad(phi, theta);

n = 0;

while norm(curr_grad)>1e-2 && n<200
    
    % Gradient ascent with step parameter 0.1
    phi = phi + 0.1*curr_grad(1);
    theta = theta + 0.1*curr_grad(2);
    
    curr_grad = grad(phi, theta);
    
    n = n + 1;
    
end

guess = polar_axis(phi, theta);

end