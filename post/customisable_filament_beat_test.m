clear;

%% Parameters

% Physical parameters
wE = 0.5*pi; % Length of the effective stroke in phi-space. wE \in (0, 2*pi)
wR = 0.5*pi; % Width of non-zero curvature region in phi-space. wR \in (0, 2*pi - wE)
wZ = 0.1*pi; % Maximum shift to phi to ensure we never have a stationary filament. Ideally should be small.
L = 1; % Filament length
A = 1.99*L; % Stroke amplitude. A \in [0, 2L]

% Numerical parameters
N = 20; % Number of points along the filament
num_phi = 100; % Number of different points in the cycle to consider

% Derived parameters
dl = L/(N-1);
beat_change_angle = acos(0.5*A/L);
dphi = 2*pi/num_phi;

%% Main loop

phi_array = dphi*(0:num_phi-1);

t_hat = zeros(2,N);
X = zeros(2*num_phi,N);

vid_fig = figure;

for m = 1:num_phi
    
    p = phi_array(m);
    
    theta = beat_change_angle + (pi - 2*beat_change_angle)*Delta(0,p,wR,wE,L)/pi;
    
    t_hat(1,1) = cos(theta);
    t_hat(2,1) = sin(theta);
    
    for n = 2:N
        
        s = (n-1)*dl;
        
        p = mod(phi_array(m) - wZ*s/L, 2*pi);
        
        theta = beat_change_angle + (pi - 2*beat_change_angle)*Delta(s,p,wR,wE,L)/pi;
        
        t_hat(1,n) = cos(theta);
        t_hat(2,n) = sin(theta);
        
        X(2*m - 1 : 2*m, n) = X(2*m - 1 : 2*m, n-1) + 0.5*dl*(t_hat(:,n-1) + t_hat(:,n));
        
    end
    
    figure(vid_fig);
    plot(X(2*m - 1,:), X(2*m,:), 'k-');
    axis equal;
    axis([-1 1 0 1]);
    pause(0.1);
    
end

figure;
hold on;
my_colours = gray(1.1*num_phi);
for m = 1:num_phi
    plot(X(2*m - 1,:), X(2*m,:), '-', 'Color', my_colours(end-m+1,:));  
end
axis equal;
axis([-1 1 0 1]);
hold off;

%% Functions

function out = Delta(s,p,w, wE,L)

phi1 = wE + s*(2*pi - wE - w)/L;

if p < wE
    
    temp = pi*(wE - p)/wE;
    
    out = temp - sin(temp)*cos(temp);
    
elseif p < phi1
    
    out = 0;
    
elseif p < phi1+w
    
    temp = pi*(p - phi1)/w;
    
    out = temp - sin(temp)*cos(temp);
    
else
    
    out = pi;
    
end

end