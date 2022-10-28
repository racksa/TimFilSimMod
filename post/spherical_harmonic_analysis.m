function spherical_harmonic_analysis(sim_name, N, frame_type)
% N is the degree at which we truncate the spherical harmonic expansion.

close all;

num_largest_modes = 10; % How many of the largest scale modes to examine.

valid_type = strcmp(frame_type, 'fixed'); % Fixed-orientation body-centred frame.
valid_type = valid_type || strcmp(frame_type, 'body'); % Body-centred frame such that a given filament base is always at the same spherical polar coordinates.
valid_type = valid_type || strcmp(frame_type, 'omega'); % Time-dependent body frame with polar axis determined by the body's angular velocity.

assert(valid_type, 'Unsupported frame type.');

%% Open files etc.

blob_references = load(strcat(sim_name, '_blob_references.dat'));
fil_references = load(strcat(sim_name, '_fil_references.dat'));
par = load(strcat(sim_name, '.par'));

NFIL = par(1);
NSEG = par(2);
NBLOB = par(4);
MU = par(5);
RSEG = par(6);
RBLOB = par(7);
DT = par(12);

L = 2.2*NSEG*RSEG;

Rbody = str2double(sim_name(10))*L;

body_state_fid = fopen(strcat(sim_name, '_body_states.dat'));
seg_state_fid = fopen(strcat(sim_name, '_seg_states.dat'));
body_vel_fid = fopen(strcat(sim_name, '_body_vels.dat'));
seg_vel_fid = fopen(strcat(sim_name, '_seg_vels.dat'));

% Assume NBOD = 1 for now
body_state_format = repmat('%f', [1 8]);
seg_state_format = repmat('%f', [1 (1 + 4*NFIL*NSEG)]);
body_vel_format = repmat('%f', [1 7]);
seg_vel_format = repmat('%f', [1 (1 + 6*NFIL*NSEG)]);

%% Form the angle grid

phi = pi*(0 : (2*N - 1))/N; % Longitudinal spacing is equidistant.

% Latitudinal spacing is given by Gauss-Legendre quadrature. This means
% that the grid points must be solved for as roots of the Legendre
% polynomial of order N+1. Since the Legendre polynomials are orthogonal,
% we can do this using the Jacobi matrix corresponding to their three-term
% recurrence relation (cf. "Calculation of Gauss Quadrature Rules" - Golub
% and Welsch [1968]).
J = zeros(N+1, N+1);

for n=1:N
    
    J(n,n+1) = n/sqrt(4*n*n - 1);
    J(n+1,n) = J(n,n+1);
    
end

[V, D] = eig(J);

theta = acos(diag(D));
w = 2*(V(1,:)).^2;

% We will also need the associated Legendre polynomials evaluated at these
% theta values.
P = zeros(N+1, N+1, N+1);

for n=0:N
    
    P(1:n+1, :, n+1) = legendre(n, cos(theta), 'norm');
    % * Each slice corresponds to a different degree.
    % * Each column corresponds to a choice of theta.
    % * Each row corresponds to a different order, which is at most the
    % degree.
    
end

%% Find the points at which to estimate the flow field

r = Rbody + 1.1*L;

flow_points = zeros(3, N+1, 2*N);
theta_dir = zeros(3, N+1, 2*N);

for i=1:N+1
    for j=1:2*N
        
        flow_points(:, i, j) = r*[cos(phi(j))*sin(theta(i)); sin(phi(j))*sin(theta(i)); cos(theta(i))];
        
        theta_dir(:, i, j) = [cos(phi(j))*cos(theta(i)); sin(phi(j))*cos(theta(i)); -sin(theta(i))];
        
    end
end

%% Form variables required for identifying the largest modes

num_largest_modes = min([num_largest_modes, 3*(N+1)*(N+2)]);

plotting_periods = 50;

scale_check_frequency = 20; % Check mode magnitudes every scale_check_frequency lines.
periods_per_scale_check = ceil(scale_check_frequency/floor(floor(54.4389/DT)/par(13)));

steady_state_periods = 300; % After how many periods do we assume the system is at steady-state behaviour.

radial_mode_scales = zeros(1, 2*(N+1)*(N+1));
radial_mode_names = cell(1, num_largest_modes);
radial_plot_lines = NaN(num_largest_modes+1, plotting_periods*ceil(ceil(54.4389/DT)/par(13)));
radial_mode_scales_plot = figure;
radial_largest_modes_plot = figure;

theta_mode_scales = zeros(1, 2*(N+1)*(N+1));
theta_mode_names = cell(1, num_largest_modes);
theta_plot_lines = NaN(num_largest_modes+1, plotting_periods*ceil(ceil(54.4389/DT)/par(13)));
theta_mode_scales_plot = figure;
theta_largest_modes_plot = figure;

phi_mode_scales = zeros(1, 2*(N+1)*(N+1));
phi_mode_names = cell(1, num_largest_modes);
phi_plot_lines = NaN(num_largest_modes+1, plotting_periods*ceil(ceil(54.4389/DT)/par(13)));
phi_mode_scales_plot = figure;
phi_largest_modes_plot = figure;

plot_pos = 1;

%% Process the data

lines_read = 0;

while ~(feof(body_state_fid) || feof(seg_state_fid) || feof(body_vel_fid) || feof(seg_vel_fid))
    
    try

        body_state_line = textscan(body_state_fid, body_state_format, 1, 'CommentStyle', '%', 'Delimiter', ' ');
        seg_state_line = textscan(seg_state_fid, seg_state_format, 1, 'CommentStyle', '%', 'Delimiter', ' ');
        body_vel_line = textscan(body_vel_fid, body_vel_format, 1, 'CommentStyle', '%', 'Delimiter', ' ');
        seg_vel_line = textscan(seg_vel_fid, seg_vel_format, 1, 'CommentStyle', '%', 'Delimiter', ' ');
        
        nt = body_state_line{1};
        T = nt*DT/54.4389;
        
        if T>=steady_state_periods
            
            q = cell2mat(body_state_line(5:8));
            qsq = q.^2;
            Q = [1 - 2*(qsq(3) + qsq(4)), 2*(q(2)*q(3) - q(1)*q(4)), 2*(q(2)*q(4) + q(1)*q(3));...
                2*(q(2)*q(3) + q(1)*q(4)), 1 - 2*(qsq(2) + qsq(4)), 2*(q(3)*q(4) - q(1)*q(2));...
                2*(q(2)*q(4) - q(1)*q(3)), 2*(q(3)*q(4) + q(1)*q(2)), 1 - 2*(qsq(2) + qsq(3))];
            % Q maps from the reference configuration to the current one.
            
            % R maps from the reference configuration to the chosen one.
            if strcmp(frame_type, 'fixed')
                
                R = Q;
                
            elseif strcmp(frame_type, 'body')
                
                R = eye(3);
                
            elseif strcmp(frame_type, 'omega')
                
                axis = Q' * cell2mat(body_vel_line(5:7))';
                axis = axis/norm(axis);
                rot_axis = cross(axis, [0;0;1]);
                R = [0, -rot_axis(3), rot_axis(2); rot_axis(3), 0, -rot_axis(1); -rot_axis(2), rot_axis(1), 0];
                R = eye(3) + R + R*R/(1 + axis(3));
                
            end
            
            U = R * Q' * cell2mat(body_vel_line(2:4))';
            Omega = R * Q' * cell2mat(body_vel_line(5:7))';
            
            %% Evaluate the approx. flow field
            
            v = zeros(3, N+1, 2*N);
            
            for n=1:NBLOB
                
                blob_pos = R * blob_references(:,n);
                
                blobU = U + cross(Omega, blob_pos);
                
                for i=1:N+1
                    for j=1:2*N
                        
                        v(:, i, j) = v(:, i, j) + flow_field(flow_points(:, i, j), blob_pos, blobU, [0;0;0], RBLOB);
                        
                    end
                end
                
            end
            
            for n=1:NFIL
                
                seg_pos = R * fil_references(:,n);
                
                for m=1:NSEG
                    
                    p = 2 + 6*((n-1)*NSEG + m - 1);
                    
                    segU = R * Q' * cell2mat(seg_vel_line(p:p+2))';
                    segOmega = R * Q' * cell2mat(seg_vel_line(p+3:p+5))';
                    
                    for i=1:N+1
                        for j=1:2*N
                            
                            v(:, i, j) = v(:, i, j) + flow_field(flow_points(:, i, j), seg_pos, segU, segOmega, RSEG);
                            
                        end
                    end
                    
                    if m<NSEG
                        
                        id = 2 + 4*((n-1)*NSEG + m - 1);
                        
                        q = cell2mat(seg_state_line(id:id+3));
                        qsq = q.^2;
                        t = [1 - 2*(qsq(3) + qsq(4)); 2*(q(2)*q(3) + q(1)*q(4)); 2*(q(2)*q(4) - q(1)*q(3))];
                        
                        id = id + 4;
                        
                        q = cell2mat(seg_state_line(id:id+3));
                        qsq = q.^2;
                        t = t + [1 - 2*(qsq(3) + qsq(4)); 2*(q(2)*q(3) + q(1)*q(4)); 2*(q(2)*q(4) - q(1)*q(3))];
                        
                        t = 1.1*RSEG*t;
                        
                        seg_pos = seg_pos + (R * Q' * t);
                        
                    end
                    
                end
                
            end
            
            %% Convert the approx. flow field to polar spherical coords
            
            v_r = zeros(N+1, 2*N);
            v_theta = zeros(N+1, 2*N);
            v_phi = zeros(N+1, 2*N);
            
            for i=1:N+1
                for j=1:2*N
                    
                    radial_dir = flow_points(:, i, j)/norm(flow_points(:, i, j));
                    v_r(i, j) = v(:, i, j)' * radial_dir;
                    
                    v_theta(i, j) = v(:, i, j)' * theta_dir(:, i, j);
                    
                    phi_dir = cross(radial_dir, theta_dir(:, i, j));
                    v_phi(i, j) = v(:, i, j)' * phi_dir;
                    
                end
            end
            
            %% Solve for the coefficients
            
            Y = fft(v_r')/N;
            
            Am = real(Y(1:N+1, :))';
            Am(:,1) = 0.5*Am(:,1);
            Am(:,N+1) = 0.5*Am(:,N+1);
            Am = (w.*Am')'; % Each row corrsponds to a different theta and each column to a different order m.
            
            Bm = -imag(Y(1:N+1, :))';
            Bm = (w.*Bm')';
            
            for n=N:-1:0
                for m=n:-1:0
                    
                    Anm(n+1,m+1) = P(m+1, :, n+1)*Am(:, m+1);
                    Bnm(n+1,m+1) = P(m+1, :, n+1)*Bm(:, m+1);
                    
                end
            end
            
            Y = fft(v_theta')/N;
            
            Cm = real(Y(1:N+1, :))';
            Cm(:,1) = 0.5*Cm(:,1);
            Cm(:,N+1) = 0.5*Cm(:,N+1);
            Cm = (w.*Cm')';
            
            Dm = -imag(Y(1:N+1, :))';
            Dm = (w.*Dm')';
            
            for n=N:-1:0
                for m=n:-1:0
                    
                    Cnm(n+1,m+1) = P(m+1, :, n+1)*Cm(:, m+1);
                    Dnm(n+1,m+1) = P(m+1, :, n+1)*Dm(:, m+1);
                    
                end
            end
            
            Y = fft(v_phi')/N;
            
            Em = real(Y(1:N+1, :))';
            Em(:,1) = 0.5*Em(:,1);
            Em(:,N+1) = 0.5*Em(:,N+1);
            Em = (w.*Em')';
            
            Fm = -imag(Y(1:N+1, :))';
            Fm = (w.*Fm')';
            
            for n=N:-1:0
                for m=n:-1:0
                    
                    Enm(n+1,m+1) = P(m+1, :, n+1)*Em(:, m+1);
                    Fnm(n+1,m+1) = P(m+1, :, n+1)*Fm(:, m+1);
                    
                end
            end
            
            %% Plot etc.
            
            radial_mode_vec = [Anm(1:end) Bnm(1:end)];
            radial_mode_scales = max([radial_mode_scales; abs(radial_mode_vec)]);
            
            theta_mode_vec = [Cnm(1:end) Dnm(1:end)];
            theta_mode_scales = max([theta_mode_scales; abs(theta_mode_vec)]);
            
            phi_mode_vec = [Enm(1:end) Fnm(1:end)];
            phi_mode_scales = max([phi_mode_scales; abs(phi_mode_vec)]);
            
            lines_read = lines_read + 1;
            
            if mod(lines_read,20)==0
                
                % v_r
                
                [~, largest_modes] = maxk(radial_mode_scales, num_largest_modes);
                
                coeff_groups = {'A', 'B'};
                
                for i=1:num_largest_modes
                    
                    group = ceil(largest_modes(i)/((N+1)*(N+1)));
                    group_mode = largest_modes(i) - (group-1)*(N+1)*(N+1);
                    m = ceil(group_mode/(N+1));
                    n = group_mode - (m-1)*(N+1);
                    
                    radial_mode_names{i} = [coeff_groups{group} '_' num2str(n-1) '^' num2str(m-1)];
                    
                end
                
                plot_mode_names = categorical(radial_mode_names);
                plot_mode_names = reordercats(plot_mode_names, radial_mode_names);
                
                figure(radial_mode_scales_plot);
                bar(plot_mode_names, radial_mode_scales(largest_modes)/radial_mode_scales(largest_modes(1)));
                title({'$v_r$', sprintf('$t/T \\in [%g, %g]$', steady_state_periods, T)}, 'Interpreter', 'latex');
                ylabel('Relative magnitude', 'Interpreter', 'latex');
                set(gca, 'FontName', 'Times', 'FontSize', 14);
                
                % v_theta
                
                [~, largest_modes] = maxk(theta_mode_scales, num_largest_modes);
                
                coeff_groups = {'C', 'D'};
                
                for i=1:num_largest_modes
                    
                    group = ceil(largest_modes(i)/((N+1)*(N+1)));
                    group_mode = largest_modes(i) - (group-1)*(N+1)*(N+1);
                    m = ceil(group_mode/(N+1));
                    n = group_mode - (m-1)*(N+1);
                    
                    theta_mode_names{i} = [coeff_groups{group} '_' num2str(n-1) '^' num2str(m-1)];
                    
                end
                
                plot_mode_names = categorical(theta_mode_names);
                plot_mode_names = reordercats(plot_mode_names, theta_mode_names);
                
                figure(theta_mode_scales_plot);
                bar(plot_mode_names, theta_mode_scales(largest_modes)/theta_mode_scales(largest_modes(1)));
                title({'$v_\theta$', sprintf('$t/T \\in [%g, %g]$', steady_state_periods, T)}, 'Interpreter', 'latex');
                ylabel('Relative magnitude', 'Interpreter', 'latex');
                set(gca, 'FontName', 'Times', 'FontSize', 14);
                
                % v_phi
                
                [~, largest_modes] = maxk(phi_mode_scales, num_largest_modes);
                
                coeff_groups = {'E', 'F'};
                
                for i=1:num_largest_modes
                    
                    group = ceil(largest_modes(i)/((N+1)*(N+1)));
                    group_mode = largest_modes(i) - (group-1)*(N+1)*(N+1);
                    m = ceil(group_mode/(N+1));
                    n = group_mode - (m-1)*(N+1);
                    
                    phi_mode_names{i} = [coeff_groups{group} '_' num2str(n-1) '^' num2str(m-1)];
                    
                end
                
                plot_mode_names = categorical(phi_mode_names);
                plot_mode_names = reordercats(plot_mode_names, phi_mode_names);
                
                figure(phi_mode_scales_plot);
                bar(plot_mode_names, phi_mode_scales(largest_modes)/phi_mode_scales(largest_modes(1)));
                title({'$v_\phi$', sprintf('$t/T \\in [%g, %g]$', steady_state_periods, T)}, 'Interpreter', 'latex');
                ylabel('Relative magnitude', 'Interpreter', 'latex');
                set(gca, 'FontName', 'Times', 'FontSize', 14);
                
                pause(0.00001);
                
            end
            
            if (T > steady_state_periods+periods_per_scale_check)&&(T <= steady_state_periods+periods_per_scale_check+plotting_periods)
                
                % v_r
                
                radial_plot_lines(1, plot_pos) = T;
                
                [~, largest_modes] = maxk(radial_mode_scales, round(num_largest_modes/2));
                
                figure(radial_largest_modes_plot);
                
                for i=1:round(num_largest_modes/2)
                    
                    radial_plot_lines(1+i, plot_pos) = radial_mode_vec(largest_modes(i));
                    
                    legend_string = ['$' radial_mode_names{i} '$'];
                    
                    plot(radial_plot_lines(1,:), radial_plot_lines(i+1,:), 'DisplayName', legend_string);
                    
                    hold on;
                    
                end
                
                hold off;
                
                legend('show', 'Interpreter', 'latex');
                xlabel('$t/T$', 'Interpreter', 'latex');
                title('$v_r$', 'Interpreter', 'latex');
                set(gca, 'FontName', 'Times', 'FontSize', 14);
                
                % v_theta
                
                theta_plot_lines(1, plot_pos) = T;
                
                [~, largest_modes] = maxk(theta_mode_scales, round(num_largest_modes/2));
                
                figure(theta_largest_modes_plot);
                
                for i=1:round(num_largest_modes/2)
                    
                    theta_plot_lines(1+i, plot_pos) = theta_mode_vec(largest_modes(i));
                    
                    legend_string = ['$' theta_mode_names{i} '$'];
                    
                    plot(theta_plot_lines(1,:), theta_plot_lines(i+1,:), 'DisplayName', legend_string);
                    
                    hold on;
                    
                end
                
                hold off;
                
                legend('show', 'Interpreter', 'latex');
                xlabel('$t/T$', 'Interpreter', 'latex');
                title('$v_\theta$', 'Interpreter', 'latex');
                set(gca, 'FontName', 'Times', 'FontSize', 14);
                
                % v_phi
                
                phi_plot_lines(1, plot_pos) = T;
                
                [~, largest_modes] = maxk(phi_mode_scales, round(num_largest_modes/2));
                
                figure(phi_largest_modes_plot);
                
                for i=1:round(num_largest_modes/2)
                    
                    phi_plot_lines(1+i, plot_pos) = phi_mode_vec(largest_modes(i));
                    
                    legend_string = ['$' phi_mode_names{i} '$'];
                    
                    plot(phi_plot_lines(1,:), phi_plot_lines(i+1,:), 'DisplayName', legend_string);
                    
                    hold on;
                    
                end
                
                hold off;
                
                legend('show', 'Interpreter', 'latex');
                xlabel('$t/T$', 'Interpreter', 'latex');
                title('$v_\phi$', 'Interpreter', 'latex');
                set(gca, 'FontName', 'Times', 'FontSize', 14);
                
                pause(0.00001);
                
                plot_pos = plot_pos + 1;
                
            end
            
            if (T > steady_state_periods+periods_per_scale_check+plotting_periods)
                
                fclose('all');
                break;
                
            end
            
        end
        
    catch error_message
        
        fprintf(['Function returned on the error "' error_message.message '"\n']);
        fclose('all');
        break;
        
    end
    
end

end