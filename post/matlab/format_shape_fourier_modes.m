% This script takes a .fourier_modes file describing the shape of an
% axisymmetric body in the same form as those accepted by the simulation
% code. It replaces this file with one of the same name, now storing the
% modes required to describe the same shape but centred at the origin and
% with unit length parallel to axis of the body about which the surface of
% revolution is formed.

clear;

shape_name = 'sphere';

D = load([shape_name '.fourier_modes']);

% Begin by scaling the existing shape to have unit length.
curr_length = radius(0, D(2:end)) + radius(pi, D(2:end));
D(2:end) = D(2:end)/curr_length;

% Identify the current centre of the shape.
[centre_estimate, cross_section_area, surface_area, volume] = find_centre(D(2:end));

% Shift the centre to the origin.
num_theta = max([80, 2*(D(1) - 1)]); % Guaranteeing that num_theta is even ensures we don't have to adjust the length again at the end.
theta = linspace(0, 2*pi*(1 - 1/num_theta), num_theta);
equiv_old_theta = zeros(size(theta));
for n = 2:num_theta % Don't bother with theta = 0 because it will always map to 0.
    if theta(n)<pi
        upper = pi;
        lower = equiv_old_theta(n-1);
        equiv_old_theta(n) = bisection_solve(theta(n), upper, lower, centre_estimate, D(2:end));
    elseif theta(n)==pi
        equiv_old_theta(n) = pi;
    else
        upper = 2*pi;
        lower = equiv_old_theta(n-1);
        equiv_old_theta(n) = bisection_solve(theta(n), upper, lower, centre_estimate, D(2:end));
    end
end

X = radius(equiv_old_theta, D(2:end)) .* [cos(equiv_old_theta); sin(equiv_old_theta)];
X = X - [centre_estimate; 0];
new_radii = sqrt(sum(X.^2, 1));

max_num_wavenumbers = 1 + ceil(0.5*(num_theta-1));
rhat = fft(new_radii)/num_theta;
cos_coeffs = 2*real(rhat(1:max_num_wavenumbers));
cos_coeffs(1) = 0.5*cos_coeffs(1);
if (max_num_wavenumbers ~= 1)&&(max_num_wavenumbers ~= 1 + 0.5*(num_theta-1))
    cos_coeffs(max_num_wavenumbers) = 0.5*cos_coeffs(max_num_wavenumbers);
end
sin_coeffs = -2*imag(rhat(1:max_num_wavenumbers));
Dnew = zeros(1, 1 + 2*max_num_wavenumbers);
Dnew(1) = max_num_wavenumbers;
Dnew(2:2:end) = cos_coeffs;
Dnew(3:2:end) = sin_coeffs;

% Write the new modes back to the file.
dlmwrite([shape_name '.fourier_modes'], Dnew, ' ');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Local function(s)
function r = radius(theta, modes)

    num_modes = size(modes,2)/2;
    cos_mat = cos((0:num_modes-1)' .* theta);
    sin_mat = sin((0:num_modes-1)' .* theta);
    r = modes(1:2:end)*cos_mat + modes(2:2:end)*sin_mat;

end

function [c, cs_area, s_area, vol] = find_centre(modes)

    num_intervals = 1000;

    theta = linspace(pi, 0, num_intervals+1);
    r = radius(theta, modes);
    x = r .* cos(theta);
    dx = x(2:end) - x(1:end-1);
    y = r .* sin(theta);
    c = x .* y;
    c = sum(dx .* (c(1:end-1) + c(2:end))); % Trapezium rule.
    cs_area = sum(dx .* (y(1:end-1) + y(2:end))); % Trapezium rule. This is the area of a cross-section containing the axis of rotation.
    c = c/cs_area;

    dydx = (y(2:end) - y(1:end-1))./dx;
    dydx = [dydx(1) 0.5*(dydx(1:end-1) + dydx(2:end)) dydx(end)];
    s_area_integrand = 2*pi*y.*sqrt(1 + dydx.^2);
    s_area = 0.5*sum(dx .* (s_area_integrand(1:end-1) + s_area_integrand(2:end))); % Trapezium rule.
    
    disc_area = pi*y.*y;
    vol = 0.5*sum(dx .* (disc_area(1:end-1) + disc_area(2:end))); % Trapezium rule.

end

function e = equiv_angle_error(new, old, centre, modes)

    r = radius(old, modes);
    e = new - mod(atan2(r*sin(old), r*cos(old) - centre), 2*pi);

end

function old = bisection_solve(new, upper, lower, centre, modes)
    
    while upper-lower > 1e-4
        
        old = 0.5*(upper + lower);
        
        error = equiv_angle_error(new, old, centre, modes);
        
        if error > 0 % old angle guess isn't big enough
            
            lower = old;
            
        else
            
            upper = old;
            
        end
        
    end
    
    old = 0.5*(upper + lower);

end