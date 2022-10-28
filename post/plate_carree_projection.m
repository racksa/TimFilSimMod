function plate_carree_projection(v3, seg_pos, NFIL, NSEG, num_shown_frames, line_width)
% Un-wrap a spherical surface using the basic "plate carr\'{e}e" equirectangular projection.
% Takes inputs from during the execution of animate.m.

% Run v3 = find_the_pole(seg_pos, NFIL, NSEG); for fixed spheres whose
% coordinated state is clearly azimuthal, polar or a linear combination of
% the two.

if abs(v3(3))==1
    
    v1 = [1;0;0];
    
else
    
    v1 = [-v3(2); v3(1); 0];
    v1 = v1/norm(v1);
    
end

v2 = cross(v3, v1);

phi = zeros(1,NFIL);
theta = zeros(1,NFIL);

R = norm(seg_pos(:,1,1));

for n=1:NFIL
    p = n + (n-1)*NSEG;
    a1 = v1' * seg_pos(:,p,end);
    a2 = v2' * seg_pos(:,p,end);
    a3 = v3' * seg_pos(:,p,end);
    phi(n) = atan2(a2,a1);
    theta(n) = acos(a3/R);
end

x_base = R*phi;
y_base = R*(0.5*pi - theta);
planar_seg_pos = NaN(size(seg_pos));

for n=1:NFIL
    
    Q = [cos(phi(n))*v2 - sin(phi(n))*v1, ...
        sin(theta(n))*v3 - cos(theta(n))*(cos(phi(n))*v1 + sin(phi(n))*v2), ...
        sin(theta(n))*(cos(phi(n))*v1 + sin(phi(n))*v2) + cos(theta(n))*v3]';
    
    for m=1:num_shown_frames
        
        p = n + (n-1)*NSEG;
        planar_seg_pos(:,p,m) = [x_base(n); y_base(n); 0];
        planar_seg_pos(:, p+1:p+NSEG-1, m) = Q*(seg_pos(:, p+1:p+NSEG-1, m) - seg_pos(:,p,m)) + planar_seg_pos(:,p,m);
        
    end
    
end

figure;
hold on;

for m=num_shown_frames:-1:1
    
    plot3(planar_seg_pos(1,:,m), planar_seg_pos(2,:,m), planar_seg_pos(3,:,m), '-', 'LineWidth', line_width, 'Color', [1 1 1]*(1 - m/num_shown_frames));
    
end

axis equal;
axis off;

% The second part of this function plots the aplanarity of the filament
% motions against their theta position, to quantify how the motion varies
% as we approach the defects.
%
% By replacing planar_seg_pos with seg_pos, this part of the function can
% also be used on planar simulations.

aplanarity = zeros(1, NFIL);
beat_angle = zeros(1, NFIL);

for n=1:NFIL
    
    p = n + NSEG*(n-1);
    
    proj = reshape(planar_seg_pos(1:2, p:p+NSEG-1, :), 2, num_shown_frames*NSEG); % project onto plane
    proj = proj - mean(proj, 2); % centre
    [~, v, pw, ~] = PCA(proj'); % pass transpose because my PCA function expects each row to be a data point
    
    aplanarity(n) = (max(v(2,:)*proj) - min(v(2,:)*proj))/(max(v(1,:)*proj) - min(v(1,:)*proj)); % bounding rectangle
%     aplanarity(n) = sqrt(pw(2)/pw(1)); % best-fit ellipse

    beat_angle(n) = atan2(v(1,2)*sign(v(1,2)), v(1,1)*sign(v(1,1))); % We want the angle between 0 and pi/2.
    
end

figure;
scatter(theta, aplanarity, 'ko');
box on;
set(gca, 'XLim', [0 pi]);
set(gca, 'FontName', 'Times', 'FontSize', 24);
xticks([0 pi/4 pi/2 3*pi/4 pi]);
xticklabels({'0', '\pi/4', '\pi/2', '3\pi/4', '\pi'});
xlabel('$\theta$', 'Interpreter', 'latex');
ylabel('$r_2/r_1$', 'Interpreter', 'latex');

figure;
scatter(theta, beat_angle, 'ko');
box on;
axis([0 pi 0 pi/2]);
set(gca, 'FontName', 'Times', 'FontSize', 24);
xticks([0 pi/4 pi/2 3*pi/4 pi]);
xticklabels({'0', '\pi/4', '\pi/2', '3\pi/4', '\pi'});
yticks([0 pi/8 pi/4 3*pi/8 pi/2]);
yticklabels({'0', '\pi/8', '\pi/4', '3\pi/8', '\pi/2'});
xlabel('$\theta$', 'Interpreter', 'latex');
ylabel('Angle between beat dir. and $\hat{e}_\phi$', 'Interpreter', 'latex');
title(sprintf('Mean angle = %g', mean(beat_angle)));

end