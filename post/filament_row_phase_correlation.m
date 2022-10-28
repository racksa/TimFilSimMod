clear;

%% Inputs
sim_name = '15_filament_test';
periods_to_use = 15;
use_actual_phase_variable = false; % Requires reading through ~60 times more data if false.

%% Load data
par = load([sim_name '.par']);
saves_per_period = round(par(14)/par(13));

if use_actual_phase_variable
    
    Ptotal = load([sim_name '_filament_phases.dat']);

else
    
    % Use tip displacement as a stand-in for the actual phase
    fid = fopen([sim_name '_seg_states.dat']);
    read_line = @(fid, line_length) textscan(fid, repmat('%f', [1 line_length]), 1, 'CommentStyle', '%', 'Delimiter', ' ');
    Ptotal = zeros(ceil(par(11)/par(13)), 1 + par(1));
    for n=1:size(Ptotal,1)
        D = cell2mat(read_line(fid, 1 + 3*par(1)*par(2)));
        Ptotal(n,1) = D(1);
        for m=1:par(1)
            Ptotal(n, m+1) = D(3*m*par(2)) - D(3 + 3*par(2)*(m-1)); % Beat is in the y-direction.
        end
    end
    Pavg = mean(Ptotal(:,2:end), 'all');
    Ptotal(:,2:end) = Ptotal(:,2:end) - Pavg;
    
end

%% Evaluate the correlation function
P = Ptotal(end+1-(periods_to_use*saves_per_period):end, 2:end);

if use_actual_phase_variable
    
    corr_fun = @(p1,p2) (1 - abs(exp(sqrt(-1)*p1) - exp(sqrt(-1)*p2)));

else
    
    corr_fun = @(p1,p2) (p1 * p2);
    
end

C = zeros(size(P));
Ccount = zeros(size(P));

for i_n = 1:size(P,1)
    for j_n = 1:size(P,2)
        
        p1 = P(i_n, j_n);
        
        for i_m = i_n:size(P,1)
            for j_m = j_n:size(P,2)
                
                p2 = P(i_m, j_m);
                
                delta_t = 1 + i_m - i_n;
                delta_x = 1 + j_m - j_n;
                
                C(delta_t, delta_x) = C(delta_t, delta_x) + corr_fun(p1,p2);
                Ccount(delta_t, delta_x) = Ccount(delta_t, delta_x) + 1;
                
            end
        end
        
    end
end

C = C ./ Ccount;

if ~use_actual_phase_variable
    C  = C/C(1); % Non-unit normalisation constant in this case.
end

figure;
contourf(0:size(P,2)-1, (0:size(P,1)-1)/saves_per_period, C, 20, 'EdgeColor', 'none');
colormap gray;
caxis([-1 1]);
colorbar;
xlabel('$\Delta x/d$', 'Interpreter', 'latex');
ylabel('$\Delta t/T$', 'Interpreter', 'latex');
set(gca, 'FontName', 'Times', 'FontSize', 24);

if ~use_actual_phase_variable
    
    fclose(fid);
    
end