% TO-DO: Make it do something for NBOD != 1.

clear
close all

FileName = 'R_equals_1L_sphere6_f_220';

show_body_vel = false;

set(groot,'DefaultFigureRenderer','painters');

video = VideoWriter(strcat(FileName,'_vel_animation2'),'MPEG-4');
open(video);

body_state_fid = fopen(strcat(FileName,'_body_states.dat'));
parameters = load(strcat(FileName,'.par'));
blob_references = load(strcat(FileName,'_blob_references.dat'));
fil_references = load(strcat(FileName,'_fil_references.dat'));
seg_vel_fid = fopen(strcat(FileName,'_seg_vels.dat'));
body_vel_fid = fopen(strcat(FileName,'_body_vels.dat'));

NFIL = parameters(1);
NSEG = parameters(2);
NBOD = parameters(3);
NBLOB = parameters(4);
RSEG = parameters(6);
RBLOB = parameters(7);
DT = parameters(12);

body_radius_over_L = str2double(FileName(10));

num_frames = 1 + floor((parameters(11)-1)/parameters(13));

body_state_format = repmat('%f',[1 (1 + NBOD*7)]);
seg_vel_format = repmat('%f',[1 (1 + NBOD*6*NFIL*NSEG)]);
body_vel_format = repmat('%f',[1 (1 + NBOD*6)]);

[theta, phi] = polar_angles(0,0,fil_references);

fraction_to_skip = 0.2;

for n=1:num_frames
    
    try
        
        body_state_line = textscan(body_state_fid, body_state_format, 1, 'CommentStyle', '%', 'Delimiter', ' ');
        seg_vel_line = textscan(seg_vel_fid, seg_vel_format, 1, 'CommentStyle', '%', 'Delimiter', ' ');
        body_vel_line = textscan(body_vel_fid, body_vel_format, 1, 'CommentStyle', '%', 'Delimiter', ' ');
        
        if n >= fraction_to_skip*num_frames
            
            q = cell2mat(body_state_line(5:8));
            qsq = q.^2;
            R = [1 - 2*(qsq(3) + qsq(4)), 2*(q(2)*q(3) - q(1)*q(4)), 2*(q(2)*q(4) + q(1)*q(3));...
                2*(q(2)*q(3) + q(1)*q(4)), 1 - 2*(qsq(2) + qsq(4)), 2*(q(3)*q(4) - q(1)*q(2));...
                2*(q(2)*q(4) - q(1)*q(3)), 2*(q(3)*q(4) + q(1)*q(2)), 1 - 2*(qsq(2) + qsq(3))];
            
            for j=NFIL:-1:1
                
                base_tangent(:,j) = R*fil_references(:,j)/norm(fil_references(:,j));
                
            end
            
            for j=NFIL:-1:1
                
                id = 1 + 6*(j*NSEG - 1);
                
                end_vel(:,j) = cell2mat(seg_vel_line(id:id+2))';
                
            end
            
            end_vel_polar = polar_velocities(0, 0, end_vel, base_tangent, R);
            
            quiver(theta, phi, end_vel_polar(2,:), end_vel_polar(3,:), 'k');
            hold on;
            
            omega = R' * cell2mat(body_vel_line(5:7))';
            norm(omega)
            [theta_omega, phi_omega] = polar_angles(0,0,omega);
            
            plot(theta_omega, phi_omega, 'r^');
            plot(pi - theta_omega, phi_omega + pi*((phi_omega <= 0) - (phi_omega > 0)), 'ro');
            hold off;
            
            axis equal;
            axis([0 pi -pi pi]);
            xticks([0 0.5*pi pi]);
            xticklabels({'0', '\pi/2', '\pi'});
            yticks([-pi -0.5*pi 0 0.5*pi pi]);
            yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
            set(gca,'FontSize',24,'FontName','Times');
            xlabel('$\theta$','Interpreter','latex');
            ylabel('$\phi$','Interpreter','latex');
            
            title_string = {sprintf('$f = %g, N_{fil} = %i, R/L = %g$', parameters(14), NFIL, body_radius_over_L), sprintf('$t/T = %.2f$', DT*(1 + (n-1)*parameters(13))/36.3107)};
            title(title_string, 'Interpreter', 'latex');
            
            writeVideo(video, getframe(gcf));
            
            fprintf('Written frame %i.\n', n);
            
        end
        
    catch
        
        close(video);
        fclose(body_state_fid);
        fclose(seg_vel_fid);
        fclose(body_vel_fid);
        return;
        
    end
    
end

close(video);
fclose(body_state_fid);
fclose(seg_vel_fid);
fclose(body_vel_fid);

function [theta_array, phi_array] = polar_angles(theta0, phi0, coords)

% Find the rotation mapping the given axis to the positive z axis.
rot_axis = cross([sin(theta0)*cos(phi0); sin(theta0)*sin(phi0); cos(theta0)], [0;0;1]);
axis_norm = norm(rot_axis);
rot_axis = rot_axis/(axis_norm + 10^-16);
q = [cos(0.5*axis_norm), sin(0.5*axis_norm)*rot_axis'];
qsq = q.^2;
R = [1 - 2*(qsq(3) + qsq(4)), 2*(q(2)*q(3) - q(1)*q(4)), 2*(q(2)*q(4) + q(1)*q(3));...
    2*(q(2)*q(3) + q(1)*q(4)), 1 - 2*(qsq(2) + qsq(4)), 2*(q(3)*q(4) - q(1)*q(2));...
    2*(q(2)*q(4) - q(1)*q(3)), 2*(q(3)*q(4) + q(1)*q(2)), 1 - 2*(qsq(2) + qsq(3))];

N = size(coords, 2);

for n=N:-1:1
    
    v = R*coords(:,n)/norm(coords(:,n));
    
    theta_array(n) = acos(v(3));
    
    if (theta_array(n)==0)||(theta_array(n)==pi)
        phi_array(n) = 0;
    else
        phi_array(n) = atan2(v(2), v(1));
    end
    
end

end

function [Vpolar] = polar_velocities(theta0, phi0, Veuclidean, radial_directions, R)

axis = [sin(theta0)*cos(phi0); sin(theta0)*sin(phi0); cos(theta0)];

N = size(Veuclidean, 2);

for n=N:-1:1
    
    theta_dir = -R*axis;
    theta_dir = theta_dir - (radial_directions(:,n)' * theta_dir)*radial_directions(:,n);
    theta_dir = theta_dir/norm(theta_dir);
    
    phi_dir = cross(radial_directions(:,n), theta_dir);
    
    Vpolar(:,n) = [radial_directions(:,n)' * Veuclidean(:,n); theta_dir' * Veuclidean(:,n); phi_dir' * Veuclidean(:,n)];
    
end

end

function [out] = coord_param(V, radial_directions, R)

N = size(V,2);

Q = zeros(1,N);

for i=1:N
    
    axis_dir = R' * cross(radial_directions(:,i), V(:,i));
    axis_dir = axis_dir/norm(axis_dir);
    
    [theta, phi] = polar_angles(0, 0, axis_dir);
    
    Vpolar = polar_velocities(theta, phi, V, radial_directions, R);
    
    for j=1:N
        
        Q(i) = Q(i) + dot(Vpolar(:,i)/norm(Vpolar(:,i)),Vpolar(:,j)/norm(Vpolar(:,j)));
        
    end
    
    Q(i) = Q(i)/N;
    
end

out = mean(Q);

end
