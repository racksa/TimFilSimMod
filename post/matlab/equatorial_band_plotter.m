function equatorial_band_plotter(sim_name, aplanarity, body_state_line, NSEG, DL, seg_state_line)
% This function is intended to be called on data produced by the script
% individual_dynamics.m and as such should NOT call any clear, close etc.
% commands. It should also be careful to plot on fresh figures so as not to
% overwrite anything I was already looking at. However, by making this a
% function I should be safe to reuse varible names if I want.

fr = load([sim_name '_fil_references.dat']);
R = norm(fr(:,1));

N = 40; % How many filaments shall we assume occupy the band?

[~, I] = mink(aplanarity, N);

q = cell2mat(body_state_line(end-3:end));
bodyQ = quaternion_matrix(q);

seg_pos = zeros(3, NSEG, N);

for n=1:N
    
    seg_pos(:,1,n) = bodyQ*fr(:,I(n));
    
end

[A, V, ~, ~] = PCA(squeeze(seg_pos(:,1,:))');

v1 = V(1,:)';
v2 = V(2,:)';
v3 = cross(v1,v2);

phi = zeros(1,N);
theta = zeros(1,N);

figure;
hold on;

for n=1:N
    
    phi(n) = atan2(A(n,2),A(n,1));
    theta(n) = acos(dot(v3, seg_pos(:,1,n))/R);
    
    q = cell2mat(seg_state_line(2 + (I(n)-1)*4*NSEG : 5 + (I(n)-1)*4*NSEG));
    t1 = quaternion_matrix(q)*[1;0;0];
    
    plot3(seg_pos(1,1,n), seg_pos(2,1,n), seg_pos(3,1,n), 'k.');
    
    for m=2:NSEG
        
        q = cell2mat(seg_state_line(2 + (I(n)-1)*4*NSEG + 4*(m-1) : 5 + (I(n)-1)*4*NSEG + 4*(m-1)));
        t2 = quaternion_matrix(q)*[1;0;0];
        
        seg_pos(:,m,n) = seg_pos(:,m-1,n) + 0.5*DL*(t1+t2);
        
        plot3(seg_pos(1,m,n), seg_pos(2,m,n), seg_pos(3,m,n), 'k.');
        
    end
    
end

axis equal;
hold off;

figure;
hold on;

for n=1:N
    
    x_base = R*phi(n);
    y_base = R*(0.5*pi - theta(n));
    
    plot3(x_base, y_base, 0, 'k.');
    
    Q = [cos(phi(n))*v2 - sin(phi(n))*v1, ...
        sin(theta(n))*v3 - cos(theta(n))*(cos(phi(n))*v1 + sin(phi(n))*v2), ...
        sin(theta(n))*(cos(phi(n))*v1 + sin(phi(n))*v2) + cos(theta(n))*v3]';
    
    for m=2:NSEG
        
        temp = Q*(seg_pos(:,m,n) - seg_pos(:,1,n));
        
        plot3(x_base + temp(1), y_base + temp(2), temp(3), 'k.');
        
    end
    
%     plot(x_base, temp(3), 'k.');
    
end

axis equal;
hold off;

end

