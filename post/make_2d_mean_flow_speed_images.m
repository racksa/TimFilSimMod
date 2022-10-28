clear;

sim_name = 'custom_beat_amplitude_1.9L_fixed_swimmer';

X = load([sim_name '_flow_locations.dat']);
X = X(:,2:end);
V = load([sim_name '_flow_velocities.dat']);
V = V(:,2:end);

Nflow = round(size(V,2)/3);
Ntime = size(V,1);

rfile = load([sim_name '_blob_references.dat']);
R = norm(rfile(1:3));

% This is all cheating based on knowing how I've set-up the flow field sim.
x = X(1, 1:3:3*sqrt(Nflow));
z = X(1, 3:3*sqrt(Nflow):end);

Vavg = zeros(sqrt(Nflow), sqrt(Nflow));

for n=1:Ntime
    for m=1:Nflow
        
        my_speed = norm(V(n, 3*m - 2 : 3*m));
        [~, my_x_id] = min(abs(X(1, 3*m - 2) - x));
        [~, my_z_id] = min(abs(X(1, 3*m) - z));
        
        Vavg(my_z_id, my_x_id) = Vavg(my_z_id, my_x_id) + my_speed;
        
    end
end

Vavg = Vavg/(R*Ntime);

figure;
imagesc([z(1) z(end)]/R, [x(1) x(end)]/R, Vavg);
set(gca, 'FontName', 'Times', 'FontSize', 24);
xlabel('$x/R$', 'Interpreter', 'latex');
ylabel('$z/R$', 'Interpreter', 'latex');
colormap jet;
bar_handle = colorbar;
ylabel(bar_handle, '$vT/R$', 'Interpreter', 'latex');
axis equal;
axis image;
hold on;
theta = linspace(0, 2*pi);
plot(cos(theta), sin(theta), 'k--');
savefig(gcf, [sim_name '_squares.fig'], 'compact');

figure;
contourf(x/R, z/R, Vavg, 100, 'LineColor', 'none');
set(gca, 'FontName', 'Times', 'FontSize', 24);
xlabel('$x/R$', 'Interpreter', 'latex');
ylabel('$z/R$', 'Interpreter', 'latex');
colormap jet;
bar_handle = colorbar;
ylabel(bar_handle, '$vT/R$', 'Interpreter', 'latex');
axis equal;
hold on;
plot(cos(theta), sin(theta), 'k--');
savefig(gcf, [sim_name '_contour.fig'], 'compact');
