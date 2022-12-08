clear
close all

FileName = 'test';

fid = fopen(strcat(FileName,'.dat'));
parameters = load(strcat(FileName,'.par'));

NFIL = parameters(1);
NSEG = parameters(2);
NBOD = parameters(3);
NBLOB = parameters(4);
TOTAL_TIME_STEPS = parameters(11);
DT = parameters(12);
PLOT_FREQUENCY_IN_STEPS = parameters(13);

if exist(strcat(FileName,'_velocities.dat'),'file')
    
    vel_fid = fopen(strcat(FileName,'_velocities.dat'));
    saved_vels = true;
    fprintf('Found saved velocity data.\n');
    
else
    
    saved_vels = false;
    fprintf('No saved velocity data found. Velocities will be approximated.\n');
    
end

L = 2.2*NSEG*parameters(6);

num_frames = 1 + floor((TOTAL_TIME_STEPS-1)/PLOT_FREQUENCY_IN_STEPS);

num_data_per_body = 7 + 4*NFIL*NSEG;

textscan_format = repmat('%f',[1 (1 + NBOD*num_data_per_body)]);

X = zeros(NBOD,num_frames);
Y = zeros(NBOD,num_frames);
Z = zeros(NBOD,num_frames);
dist = zeros(NBOD,num_frames);
speed = NaN(NBOD,num_frames);

% Xg = zeros(NBOD,num_frames);
% Yg = zeros(NBOD,num_frames);
% Zg = zeros(NBOD,num_frames);

t = DT*(1:PLOT_FREQUENCY_IN_STEPS:TOTAL_TIME_STEPS)/54.4389;

for i=NBOD:-1:1
    
    coord_fig(i) = figure;
    
end

dist_fig = figure;
speed_fig = figure;

for n=1:num_frames
    
    try
        
        D = cell2mat(textscan(fid, textscan_format, 1, 'CommentStyle', '%', 'Delimiter', ' '));
        
        if saved_vels
            
            Dvel = cell2mat(textscan(vel_fid, repmat('%f',[1 (1 + NBOD*6*(1 + NFIL*NSEG))]), 1, 'CommentStyle', '%', 'Delimiter', ' '));
            
        end
        
        for i=NBOD:-1:1
            
            Dbod = D(2 + (i-1)*num_data_per_body : 1 + i*num_data_per_body);
            
            if n==1
                
                initial_pos(:,i) = Dbod(1:3)';
                
            end
            
            X(i,n) = Dbod(1) - initial_pos(1,i);
            Y(i,n) = Dbod(2) - initial_pos(2,i);
            Z(i,n) = Dbod(3) - initial_pos(3,i);
            
            dist(i,n) = norm(Dbod(1:3)' - initial_pos(:,i));
            
            if saved_vels
                
                Dvel_bod = Dvel(2 + (i-1)*6*(1 + NFIL*NSEG) : 1 + i*6*(1 + NFIL*NSEG));
                
                speed(i,n) = norm(Dvel_bod(1:3));
                
            else
                
                if n>1
                    
                    speed(i,n) = norm([X(i,n) - X(i,n-1); Y(i,n) - Y(i,n-1); Z(i,n) - Z(i,n-1)])/(DT*PLOT_FREQUENCY_IN_STEPS);
                    
                end
                
            end
            
        end
        
        fprintf('Read frame %i.\n', n);
        
        if mod(n,100)==0
            
            figure(dist_fig);
            
            for i=1:NBOD
                
                plot(t(1:n-1), dist(i, 1:n-1)/L, 'DisplayName', sprintf('Body %i',i));
                hold on;
                
            end
            
            legend1 = legend(gca,'show');
            legend1.Interpreter = 'latex';
            
            set(gca,'FontSize',24,'FontName','Times');
            xlabel('$t/T$','Interpreter','latex');
            ylabel('$d/L$','Interpreter','latex');
            
            hold off;
            
            figure(speed_fig);
            
            for i=1:NBOD
                
                plot(t(1:n-1), 54.4389*speed(i, 1:n-1)/L, 'DisplayName', sprintf('Body %i',i));
                hold on;
                
            end
            
            legend1 = legend(gca,'show');
            legend1.Interpreter = 'latex';
            
            set(gca,'FontSize',24,'FontName','Times');
            xlabel('$t/T$','Interpreter','latex');
            ylabel('$VT/L$','Interpreter','latex');
            
            hold off;
            
            for i=1:NBOD
                
                figure(coord_fig(i));
                
                plot(t(1:n-1), X(i, 1:n-1)/L, 'k-', 'DisplayName', '$(x - x_0)/L$');
                hold on;
                plot(t(1:n-1), Y(i, 1:n-1)/L, 'k--', 'DisplayName', '$(y - y_0)/L$');
                plot(t(1:n-1), Z(i, 1:n-1)/L, 'k:', 'DisplayName', '$(z - z_0)/L$');
                
                legend1 = legend(gca,'show');
                legend1.Interpreter = 'latex';
                
                set(gca,'FontSize',24,'FontName','Times');
                xlabel('$t/T$','Interpreter','latex');
                title(sprintf('Body %i', i), 'Interpreter', 'latex');
                
                hold off;
                
            end
            
            pause(0.1);
            
        end
        
    catch
        
        fclose(fid);
        
        figure(dist_fig);
        
        for i=1:NBOD
            
            plot(t(1:n-1), dist(i, 1:n-1)/L, 'DisplayName', sprintf('Body %i',i));
            hold on;
            
        end
        
        legend1 = legend(gca,'show');
        legend1.Interpreter = 'latex';
        
        set(gca,'FontSize',24,'FontName','Times');
        xlabel('$t/T$','Interpreter','latex');
        ylabel('$d/L$','Interpreter','latex');
        
        hold off;
        
        savefig(dist_fig, strcat(FileName, '_swimming_dists'));
        
        figure(speed_fig);
        
        for i=1:NBOD
            
            plot(t(1:n-1), 54.4389*speed(i, 1:n-1)/L, 'DisplayName', sprintf('Body %i',i));
            hold on;
            
        end
        
        legend1 = legend(gca,'show');
        legend1.Interpreter = 'latex';
        
        set(gca,'FontSize',24,'FontName','Times');
        xlabel('$t/T$','Interpreter','latex');
        ylabel('$VT/L$','Interpreter','latex');
        
        hold off;
        
        savefig(dist_fig, strcat(FileName, '_swimming_speeds'));
        
        for i=1:NBOD
            
            figure(coord_fig(i));
            
            plot(t(1:n-1), X(i, 1:n-1)/L, 'k-', 'DisplayName', '$(x - x_0)/L$');
            hold on;
            plot(t(1:n-1), Y(i, 1:n-1)/L, 'k--', 'DisplayName', '$(y - y_0)/L$');
            plot(t(1:n-1), Z(i, 1:n-1)/L, 'k:', 'DisplayName', '$(z - z_0)/L$');
            
            legend1 = legend(gca,'show');
            legend1.Interpreter = 'latex';
            
            set(gca,'FontSize',24,'FontName','Times');
            xlabel('$t/T$','Interpreter','latex');
            title(sprintf('Body %i', i), 'Interpreter', 'latex');
            
            hold off;
            
            savefig(coord_fig(i), strcat(FileName, sprintf('_body_%i_coords', i)));
            
        end
        
        return;
        
    end
    
end

fclose(fid);

figure(dist_fig);

for i=1:NBOD
    
    plot(t, dist/L, 'DisplayName', sprintf('Body %i',i));
    hold on;
    
end

legend1 = legend(gca,'show');
legend1.Interpreter = 'latex';

set(gca,'FontSize',24,'FontName','Times');
xlabel('$t/T$','Interpreter','latex');
ylabel('$d/L$','Interpreter','latex');

hold off;

savefig(dist_fig, strcat(FileName, '_swimming_dists'));

figure(speed_fig);

for i=1:NBOD
    
    plot(t, 54.4389*speed/L, 'DisplayName', sprintf('Body %i',i));
    hold on;
    
end

legend1 = legend(gca,'show');
legend1.Interpreter = 'latex';

set(gca,'FontSize',24,'FontName','Times');
xlabel('$t/T$','Interpreter','latex');
ylabel('$VT/L$','Interpreter','latex');

hold off;

savefig(dist_fig, strcat(FileName, '_swimming_speeds'));

for i=1:NBOD
    
    figure(coord_fig(i));
    
    plot(t, X/L, 'k-', 'DisplayName', '$(x - x_0)/L$');
    hold on;
    plot(t, Y/L, 'k--', 'DisplayName', '$(y - y_0)/L$');
    plot(t, Z/L, 'k:', 'DisplayName', '$(z - z_0)/L$');
    
    legend1 = legend(gca,'show');
    legend1.Interpreter = 'latex';
    
    set(gca,'FontSize',24,'FontName','Times');
    xlabel('$t/T$','Interpreter','latex');
    title(sprintf('Body %i', i), 'Interpreter', 'latex');
    
    hold off;
    
    savefig(coord_fig(i), strcat(FileName, sprintf('_body_%i_coords', i)));
    
end
