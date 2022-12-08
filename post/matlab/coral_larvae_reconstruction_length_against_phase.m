clear;

file_name = 'trace_realunits_1';

A1 = load([file_name '_x_cos_coeffs.dat'])';

A2 = load([file_name '_y_cos_coeffs.dat'])';

B1 = load([file_name '_x_sin_coeffs.dat'])';

B2 = load([file_name '_y_sin_coeffs.dat'])';

s_vec = @(s) [1 2*s 3*s*s];
cos_vec = @(p) cos(p*(0:(size(A1,2)-1))');
sin_vec = @(p) sin(p*(0:(size(A1,2)-1))');
 
tangent_fun = @(s,p) [s_vec(s)*(A1*cos_vec(p) + B1*sin_vec(p));...
    s_vec(s)*(A2*cos_vec(p) + B2*sin_vec(p))];

length_integrand = @(s,p) norm(tangent_fun(s,p));

num_phi = 200;

phi = (0:num_phi-1)*2*pi/num_phi;

length_array = zeros(1,num_phi);

for n = 1:num_phi
    
    s_points = 200;
    
    % Trapezium rule
    length_array(n) = 0.5*(length_integrand(0, phi(n)) + length_integrand(1, phi(n)));
    
    for m = 2:s_points-1
        
        sm = (m-1)/(s_points-1);
        
        length_array(n) = length_array(n) + length_integrand(sm, phi(n));
        
    end
    
    ds = 1/(s_points - 1);
    
    length_array(n) = ds*length_array(n);
    
end

figure;
plot(phi, length_array, 'k-');
xlabel('$\phi$', 'Interpreter', 'latex');
xticks([0 pi/2 pi 3*pi/2 2*pi]);
xticklabels({'0' '\pi/2' '\pi' '3\pi/2' '2\pi'});
set(gca, 'XLim', [0 2*pi]);
ylabel('$L(\phi)/L$', 'Interpreter', 'latex');
set(gca, 'FontName', 'Times', 'FontSize', 24);