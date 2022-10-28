% close all;

sim_name = 'single_planar_filament_f_220';
start_period = 10;
num_periods = 4;

parameters = load([sim_name '.par']);

NFIL = parameters(1);
NSEG = parameters(2);
NBOD = parameters(3);
NBLOB = parameters(4);
RSEG = parameters(6);
RBLOB = parameters(7);
DT = parameters(12);
DL = 2.2*RSEG;

seg_state_fid = fopen([sim_name '_seg_states.dat']);
seg_force_fid = fopen([sim_name '_seg_forces.dat']);
seg_vel_fid = fopen([sim_name '_seg_vels.dat']);

seg_state_format = repmat('%f', [1 (1 + 4*NBOD*NFIL*NSEG)]);
seg_force_format = repmat('%f', [1 (1 + 6*NBOD*NFIL*NSEG)]);
seg_vel_format = repmat('%f', [1 (1 + 6*NBOD*NFIL*NSEG)]);

T = 36.3107;

steps_per_period = round(T/DT);
saves_per_period = ceil(steps_per_period/parameters(13));

saves_to_use = ceil(num_periods*saves_per_period);

seg_projections = zeros(saves_to_use*(NSEG-1), 2, NFIL);

seg_tangents = zeros(saves_to_use*(NSEG-1), 3, NFIL);

distal_end_pos = zeros(3, saves_to_use, NFIL);

saves_read = 0;

while ~(feof(seg_state_fid) && feof(seg_force_fid) && feof(seg_vel_fid))
    
    seg_state_line = textscan(seg_state_fid, seg_state_format, 1, 'CommentStyle', '%', 'Delimiter', ' ');
    
    t_over_T = DT*seg_state_line{1}/T;
    
    if t_over_T >= start_period
        
        if t_over_T > start_period + num_periods
            
            aplanarity = zeros(1, NFIL);
            
            period = zeros(1, NFIL);
            
            curv_A = zeros(saves_read, NSEG-1, NFIL);
            curv_V = zeros(NSEG-1, NSEG-1, NFIL);
            curv_pw = zeros(NSEG-1, NFIL);
            curv_pce = zeros(NSEG-1, NFIL);
            
            for n=1:NFIL
                
                [~, V, pw, ~] = PCA(seg_projections(1:saves_read*(NSEG-1), :, n) - mean(seg_projections(1:saves_read*(NSEG-1), :, n),1));
                
                aplanarity(n) = pw(2)/pw(1);
                
                period(n) = find_period(distal_end_pos(:, 1:saves_read, n))/saves_per_period;
                
                for r=1:saves_read
                    
                    t1= [1 0 0];
                    
                    for m=2:NSEG
                        
                        t2 = seg_tangents((r-1)*(NSEG-1) + m - 1, :, n);
                        
                        theta1 = atan2(t1(1), t1(2)*V(1,1) + t1(3)*V(1,2));
                        theta2 = atan2(t2(1), t2(2)*V(1,1) + t2(3)*V(1,2));
                        
                        if abs(theta2-theta1)>pi
                            
                            if theta2<0
                                
                                theta2 = theta2 + 2*pi;
                                
                            else
                                
                                theta1 = theta1 + 2*pi;
                                
                            end
                            
                        end
                        
                        curv_A(r, m-1, n) = (theta2 - theta1)/DL;
                        
                        t1 = t2;
                        
                    end
                    
                end
                
                [curv_A(:,:,n), curv_V(:,:,n), curv_pw(:,n), curv_pce(:,n)] = PCA(curv_A(:,:,n));
                
            end
            
            figure;
            histogram(aplanarity, 'BinWidth', 0.05);
            
            figure;
            histogram(period, 'BinWidth', 0.05);
            
            break;
            
        end
        
        saves_read = saves_read + 1;
        
        for n=1:NFIL
            
            q = cell2mat(seg_state_line(2 + (n-1)*4*NSEG : 5 + (n-1)*4*NSEG));
            
            Qbase = quaternion_matrix(q);
            
            t1 = [1;0;0];
            
            seg_pos = [0;0;0];
            
            for m=2:NSEG
                
                q = cell2mat(seg_state_line(2 + (n-1)*4*NSEG + 4*(m-1) : 5 + (n-1)*4*NSEG + 4*(m-1)));
                
                t2 = Qbase'*quaternion_matrix(q)*[1;0;0];
                
                seg_pos = seg_pos + 0.5*DL*(t1 + t2);
                
                row = (saves_read-1)*(NSEG-1) + m - 1;
                
                seg_projections(row, 1, n) = seg_pos(2);
                seg_projections(row, 2, n) = seg_pos(3);
                
                seg_tangents(row, :, n) = t2';
                
                t1 = t2;
                
            end
            
            distal_end_pos(:, saves_read, n) = seg_pos;
            
        end
        
    end
    
end

fclose('all');