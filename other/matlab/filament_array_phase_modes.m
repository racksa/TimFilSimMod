clear;

% Provide inputs
sim_name = 'double_15_filament_test';
tol = 5; % Percentage error to allow between the filament phases and their 'low'-dimensional representation, as measured using the Frobenius norm.
t_frac = 1; % Fraction of simulation time to include, measured from the end rather than the start -- if we're omitting any data, we want it to be early, transition states.

% Load in the phases
P = load([sim_name '_filament_phases.dat']);
P = P(1 + round((1-t_frac)*(size(P,1)-1)):end,:);
t = P(:,1);
P = P(:,2:end);

% Shift the phases into [-pi, pi)
P = mod(P, 2*pi);
P(P > pi) = P(P > pi) - 2*pi;

% Perform the PCA
[A, V, pw, pce] = PCA(P);

% Identify the number of modes we consider non-negligible
for n=1:size(V,1)
    
    if pce(n)<tol
        
        num_modes_to_use = n;
        break;
        
    end
    
end

fprintf('Keeping %i/%i modes.\n', num_modes_to_use, size(V,1));

% Transfer magnitude between each mode and its coefficient.
% This ensures that each mode can be viewed as a distribution of phases
% over the array. The orthogonality of the modes is unaffected, but of
% course they're no longer orthonormal.
for n=1:num_modes_to_use
    
    [~, max_id] = max(abs(V(n,:)));
    fac = -pi/V(n,max_id);
    V(n,:) = fac*V(n,:);
    A(:,n) = A(:,n)/fac;
    
end

% Plot the modes so we can visualise them and see if they correspond to
% MCWs etc.
ref = load([sim_name '_fil_references.dat']);
x = ref(1:3:end);
y = ref(2:3:end);
[X,Y] = meshgrid(x,y);

for n=1:num_modes_to_use
    
    Z = griddata(x,y,V(n,:),X,Y);
    
    figure;
    contourf(X, Y, Z, 'LineColor', 'none');
    colormap hsv;
    caxis([-pi pi]);
    cbh = colorbar;
    cbh.Ticks = [-pi -pi/2 0 pi/2 pi];
    cbh.TickLabels = {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'};
    title(sprintf('Mode %i', n));
    axis equal;
    
end