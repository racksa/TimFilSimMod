function coeffs = az_polar_coeffs(polar_axis, seg_pos, NFIL, NSEG)
% Find the least-squares decomposition of the collection of filament
% beating directions into purely azimuthal and purely polar parts. We must
% supply the polar axis about/towards which we think the coordination is
% occuring. % Takes inputs from during the execution of animate.m.

polar_axis = polar_axis/norm(polar_axis); % Just to be sure
polar_axis = [polar_axis(1); polar_axis(2); polar_axis(3)]; % Ensure it's a column vector, not a row

if abs(polar_axis(3))==1
    
    v1 = [1;0;0];
    
else
    
    v1 = [-polar_axis(2); polar_axis(1); 0];
    v1 = v1/norm(v1);
    
end

v2 = cross(polar_axis, v1);

R = norm(seg_pos(:,1,1));

az_coeff = 0;
polar_coeff = 0;

for n=1:NFIL
    
    p = n + (n-1)*NSEG;
    
    a1 = v1' * seg_pos(:,p,end);
    a2 = v2' * seg_pos(:,p,end);
    a3 = polar_axis' * seg_pos(:,p,end);
    
    phi = atan2(a2,a1);
    theta = acos(a3/R);
    
    az_dir = cos(phi)*v2 - sin(phi)*v1;
    polar_dir = sin(theta)*polar_axis - cos(theta)*(cos(phi)*v1 + sin(phi)*v2);
    
    beat_dir = seg_pos(:,p+NSEG-1,end) - seg_pos(:,p,end);
    beat_dir = beat_dir - (beat_dir'*seg_pos(:,p,end)/norm(seg_pos(:,p,end)))*seg_pos(:,p,end)/norm(seg_pos(:,p,end)); % Remove radial component
    beat_dir = beat_dir/norm(beat_dir);
    
    az_coeff = az_coeff + az_dir'*beat_dir;
    polar_coeff = polar_coeff + polar_dir'*beat_dir;
    
end

az_coeff = az_coeff/NFIL;
polar_coeff = polar_coeff/NFIL;

coeffs = [az_coeff, polar_coeff];

end