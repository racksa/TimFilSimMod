clear
close all

FileName = '3L_over_4_spaced_wall';

fid = fopen(strcat(FileName,'.dat'));
parameters = load(strcat(FileName,'.par'));
fil_references = load(strcat(FileName,'_fil_references.dat'));

NFIL = parameters(1);
NSEG = parameters(2);
NBOD = parameters(3);
NBLOB = parameters(4);
RSEG = parameters(6);
TOTAL_TIME_STEPS = parameters(11);
DT = parameters(12);
PLOT_FREQUENCY_IN_STEPS = parameters(13);

wrong_sim_type = exist(strcat(FileName,'_blob_references.dat'),'file') || ((NBOD~=1) || (NBLOB~=0));

assert(~wrong_sim_type, 'Found parameters incompatible with RPY wall type simulation or blob reference data. This script is designed for arrays of cilia on a no-slip RPY wall.');

L = 2.2*NSEG*RSEG;

num_frames = 1 + floor((TOTAL_TIME_STEPS-1)/PLOT_FREQUENCY_IN_STEPS);

num_data_per_body = 7 + 4*NFIL*NSEG;

textscan_format = repmat('%f',[1 (1 + NBOD*num_data_per_body)]);

set(groot,'DefaultFigureRenderer','painters');

video = VideoWriter(strcat(FileName,'_wave_animation'),'MPEG-4');
open(video);

Z = zeros(1,NFIL);

time = DT*(1:PLOT_FREQUENCY_IN_STEPS:TOTAL_TIME_STEPS)/54.4389;

for n=1:num_frames
    
    try
        
        D = cell2mat(textscan(fid, textscan_format, 1, 'CommentStyle', '%', 'Delimiter', ' '));
        
        q = D(5:8);
        qsq = q.^2;
        R = [1 - 2*(qsq(3) + qsq(4)), 2*(q(2)*q(3) - q(1)*q(4)), 2*(q(2)*q(4) + q(1)*q(3));...
            2*(q(2)*q(3) + q(1)*q(4)), 1 - 2*(qsq(2) + qsq(4)), 2*(q(3)*q(4) - q(1)*q(2));...
            2*(q(2)*q(4) - q(1)*q(3)), 2*(q(3)*q(4) + q(1)*q(2)), 1 - 2*(qsq(2) + qsq(3))];
        
        
        Xc = D(2:4)';
        
        for j=NFIL:-1:1
            
            end_pos = Xc + R*fil_references(:,j);
            
            for k=2:NSEG
                
                id = 9 + 4*((j-1)*NSEG + k-2);
                
                q = D(id:id+3);
                qsq = q.^2;
                t = [1 - 2*(qsq(3) + qsq(4)); 2*(q(2)*q(3) + q(1)*q(4)); 2*(q(2)*q(4) - q(1)*q(3))];
                
                id = id + 4;
                
                q = D(id:id+3);
                qsq = q.^2;
                t = t + [1 - 2*(qsq(3) + qsq(4)); 2*(q(2)*q(3) + q(1)*q(4)); 2*(q(2)*q(4) - q(1)*q(3))];
                
                t = 1.1*RSEG*t;
                
                end_pos = end_pos + t;
                
            end
            
            Z(j) = end_pos(3);
            
        end
        
        plot(Z/L, 'k-');
        set(gca, 'Ylim', [0 1]);
        set(gca,'FontSize',24,'FontName','Times');
        xlabel('$n$','Interpreter','latex');
        ylabel('$z_N/L$','Interpreter','latex');
        title(sprintf('$t/T = %.3f$', time(n)),'Interpreter','latex');
        
        writeVideo(video, getframe(gcf));
        
        fprintf('Written frame %i.\n', n);
        
    catch
        
        close(video);
        fclose(fid);
        return;
        
    end
    
end

close(video);
fclose(fid);
