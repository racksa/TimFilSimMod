clear;

sims = dir('../bab_pair_2*filament_phases.dat');

dphi_fig = figure;
hold on;

osc_mag_fig = figure;
hold on;

for n=1:length(sims)
    
    P = load([sims(n).folder, '/', sims(n).name]);
    SIM_NAME = [sims(n).folder, '/', sims(n).name(1:end-20)];
    
    par = read_parameter_file(SIM_NAME);
    
     if isfield(par, 'TORSIONAL_SPRING_MAGNITUDE_FACTOR')
        
        legend_entry = sprintf('$k/\\overline{Q_\\theta} = %g$', par.TORSIONAL_SPRING_MAGNITUDE_FACTOR);
        
    else
        
        legend_entry = '$k/\overline{Q_\theta} = \infty$';
        par.TORSIONAL_SPRING_MAGNITUDE_FACTOR = inf;
        
     end
    
    figure(dphi_fig);
    plot(P(:,1)/par.STEPS_PER_PERIOD, P(:,3)-P(:,2), 'DisplayName', legend_entry);
    
    figure(osc_mag_fig);
    start_id = round(0.9*size(P,1));
    plot(1/par.TORSIONAL_SPRING_MAGNITUDE_FACTOR, mean(abs(P(start_id:end,3)-P(start_id:end,2))), '*', 'DisplayName', legend_entry);
    
end

box on;
xlabel('$t/T$', 'Interpreter', 'latex');
ylabel('$\Delta\varphi$', 'Interpreter', 'latex');
legend('show', 'Interpreter', 'latex');
set(gca, 'FontName', 'Times', 'FontSize', 16);