clear;

sim_name = 'MCW_parameter_random_ICs_2';

periods_to_use = 15; % Measured from the end

% Form the sequence of Kuramoto order parameters:
P = load([sim_name '_filament_phases.dat']);
par = load([sim_name '.par']);
data_per_period = round(par(14)/par(13));
start_entry = size(P,1) + 1 - periods_to_use*data_per_period;
t = P(start_entry:end,1)/par(14);
P = P(start_entry:end,2:end);
R = abs(mean(exp(sqrt(-1)*P), 2));

% Plot the return map:
figure;
plot(R(1:end-1), R(2:end), 'k-');
axis equal;
box on;
xlabel('$R_{n}$', 'Interpreter', 'latex');
ylabel('$R_{n+1}$', 'Interpreter', 'latex');
set(gca, 'FontName', 'Times', 'FontSize', 24);

% Plot the sine-of-phase contour for the same data, for comparison:
C = sin(P);
figure;
contourf(1:size(P,2), t, C, 'EdgeColor', 'none');
xlabel('$m$', 'Interpreter', 'latex');
ylabel('$t_n/T$', 'Interpreter', 'latex');
colormap 'gray';
caxis([-1 1]);
cbar_ax = colorbar;
ylabel(cbar_ax, '$\sin(\varphi)$', 'Interpreter', 'latex');
set(gca, 'FontName', 'Times', 'FontSize', 24);