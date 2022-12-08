clear

cd data;

FileName = 'planar_test_14';%'curved_plane_10_by_14_R_equals_-10L';

fid = fopen(strcat(FileName,'.dat'));
parameters = load(strcat(FileName,'.par'));

Nf = parameters(1);
NperFil = parameters(2);
Nframes = floor(parameters(9)/10);
dt = parameters(10);

textscan_format = repmat('%f',[1 (1 + 13*Nf*NperFil)]);

V = zeros(1,Nframes);
Q = NaN(1,Nframes);
Dx_instant = zeros(Nf,Nframes);
Dx = NaN(1,Nframes);
Dy_instant = zeros(Nf,Nframes);
Dy = NaN(1,Nframes);

frames_per_period = ceil(36.5976/(10*dt));

avgFrames = frames_per_period;

Vtemp = zeros(3,Nf);

Q_plot = figure;
direction_plot = figure;

for frame_num = 1:Nframes
    
    D = cell2mat(textscan(fid,textscan_format,1,'CommentStyle','%','Delimiter',' '));
    
    for i=1:Nf
        
        base_id = 2 + 13*NperFil*(i-1);
        id = 2 + 13*(i*NperFil - 1);
        
        beat_pos = D(id:id+2) - D(base_id:base_id+2);
        
        q = D(base_id+3:base_id+6);
        qsq = q.^2;
        
        d3 = [1 - 2*(qsq(3) + qsq(4)), 2*(q(2)*q(3) + q(1)*q(4)), 2*(q(2)*q(4) - q(1)*q(3))];
        
        d1 = [1 0 0];
        d2 = cross(d3, d1);
        
        Vtemp(1,i) = dot(beat_pos, d1);
        Vtemp(2,i) = dot(beat_pos, d2);
        Vtemp(3,i) = dot(beat_pos, d3);
        
        Dx_instant(i,frame_num) = abs(Vtemp(1,i))/norm(Vtemp(1:2,i));
        Dy_instant(i,frame_num) = abs(Vtemp(2,i))/norm(Vtemp(1:2,i));
        
    end
    
    V(frame_num) = norm(cov(Vtemp'),'fro');
    
    if frame_num >= avgFrames
        
        Q(frame_num) = mean(V(1+frame_num-avgFrames:frame_num));
        
        Dxn = mean(Dx_instant(:,1+frame_num-avgFrames:frame_num), 2);
        Dx(frame_num) = mean(Dxn);
        
        Dyn = mean(Dy_instant(:,1+frame_num-avgFrames:frame_num), 2);
        Dy(frame_num) = mean(Dyn);
        
        phin = atan(Dyn./Dxn);
        phi = mean(phin)
        
    end
    
    figure(Q_plot);
    plot((1:frame_num)/frames_per_period, Q(1:frame_num), 'k-');
    set(gca,'FontSize',24,'FontName','Times');
    xlabel('$t/T$','Interpreter','latex');
    ylabel('$Q(t)$','Interpreter','latex');
    
    figure(direction_plot);
    plot((1:frame_num)/frames_per_period, Dx(1:frame_num), 'k--', 'DisplayName', '$D_x(t)$');
    hold on;
    plot((1:frame_num)/frames_per_period, Dy(1:frame_num), 'k-', 'DisplayName', '$D_y(t)$');
    hold off;
    set(gca,'FontSize',24,'FontName','Times');
    xlabel('$t/T$','Interpreter','latex');
    legend1 = legend(gca,'show');
    legend1.Interpreter = 'latex';
    
    pause(0.1);
    
end

cd ..