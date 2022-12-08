clear;

A1 = [-0.327 0.393 -0.097 0.079; 0.3935 -1.516 0.032 -0.302; 0.101 0.716 -0.118 0.142];

A2 = [0.9475 -0.018 0.158 0.01; -0.276 -0.126 -0.341 0.035; 0.048 0.263 0.186 -0.067];

B1 = [0 0.284 0.006 -0.059; 0 1.045 0.317 0.226; 0 -1.017 -0.276 -0.196];

B2 = [0 0.192 -0.05 0.012; 0 -0.499 0.423 0.138; 0 0.339 -0.327 -0.114];

s_vec = @(s) [1 2*s 3*s*s];
 
tangent_fun = @(s,p) [s_vec(s)*(A1*[1; cos(p); cos(2*p); cos(3*p)] + B1*[0; sin(p); sin(2*p); sin(3*p)]); s_vec(s)*(A2*[1; cos(p); cos(2*p); cos(3*p)] + B2*[0; sin(p); sin(2*p); sin(3*p)])];

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