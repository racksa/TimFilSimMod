D = load(['isolated_filament_for_array_beat_comparison_f_', num2str(parameters(14)),'_seg_states.dat']);
j = round(NFIL/2);
p = j + (j-1)*NSEG;
array_seg_pos = seg_pos(:, p:p+NSEG-1, :);
array_seg_pos = reshape(array_seg_pos, 3, NSEG*num_shown_frames);

% Identify principal direction of filament beat shape
[~, V, ~, ~] = PCA(array_seg_pos(1:2, :)' - mean(array_seg_pos(1:2, :), 2)');
dir = [V(1,:)'; 0];

% Construct the rotation mapping dir to the x-axis
theta = atan2(dir(2), dir(1));
q = [cos(-0.5*theta), sin(-0.5*theta)*[0 0 1]]; % rotate through -theta about the z-axis
Q = quaternion_matrix(q);

% Use this to align the filament shape with the isolated filament
isolated_base = load('isolated_filament_for_array_beat_comparison_fil_references.dat');
array_seg_pos = isolated_base + Q*(array_seg_pos - array_seg_pos(:,1));

% Compare the shapes, constructing the required isolated shapes as we go
min_error = 10^10;
min_error_id = 0;
freq_in_steps = parameters(13);
error_scale = norm(2.2*RSEG*ones(size(array_seg_pos)), 'fro'); % error we would get if every segment was out by a segment length

for i=1:num_shown_frames*freq_in_steps
    
    % Construct the isolated shapes given that we start at shape i
    isolated_seg_pos = zeros(size(array_seg_pos));
    
    for j=1:num_shown_frames
        
        row_of_D = i + (j-1)*freq_in_steps;
        
        p = 1 + (j-1)*NSEG;
        
        isolated_seg_pos(:,p) = isolated_base;
        
        for k=2:NSEG
            
            id = 4*(k-2) + 2;
            
            q = D(row_of_D, id:id+3);
            qsq = q.^2;
            t = [1 - 2*(qsq(3) + qsq(4)); 2*(q(2)*q(3) + q(1)*q(4)); 2*(q(2)*q(4) - q(1)*q(3))];
            
            id = id + 4;
            
            q = D(row_of_D, id:id+3);
            qsq = q.^2;
            t = t + [1 - 2*(qsq(3) + qsq(4)); 2*(q(2)*q(3) + q(1)*q(4)); 2*(q(2)*q(4) - q(1)*q(3))];
            
            t = 1.1*RSEG*t;
            
            p = p + 1;
            
            isolated_seg_pos(:,p) = isolated_seg_pos(:,p-1) + t;
            
        end
        
    end
    
    error = norm(array_seg_pos - isolated_seg_pos, 'fro')/error_scale;
    
    if error < min_error
        
        min_error = error;
        min_error_id = i;
        
    end
    
end

% Plot the shape sequence that gave the smallest error
figure;
i = min_error_id;
isolated_seg_pos = zeros(size(array_seg_pos));
num_plot_rows = round(sqrt(num_shown_frames));
num_plot_cols = ceil(num_shown_frames/num_plot_rows);

for j=1:num_shown_frames
    
    subplot(num_plot_rows, num_plot_cols, j);
    view(0,0);
    hold on;
    
    plot3((j-1)*L + array_seg_pos(1, 1 + (j-1)*NSEG: j*NSEG), array_seg_pos(2, 1 + (j-1)*NSEG: j*NSEG), array_seg_pos(3, 1 + (j-1)*NSEG: j*NSEG), 'k-');
    
    row_of_D = i + (j-1)*freq_in_steps;
    
    p = 1 + (j-1)*NSEG;
    
    isolated_seg_pos(:,p) = isolated_base;
    
    for k=2:NSEG
        
        id = 4*(k-2) + 2;
        
        q = D(row_of_D, id:id+3);
        qsq = q.^2;
        t = [1 - 2*(qsq(3) + qsq(4)); 2*(q(2)*q(3) + q(1)*q(4)); 2*(q(2)*q(4) - q(1)*q(3))];
        
        id = id + 4;
        
        q = D(row_of_D, id:id+3);
        qsq = q.^2;
        t = t + [1 - 2*(qsq(3) + qsq(4)); 2*(q(2)*q(3) + q(1)*q(4)); 2*(q(2)*q(4) - q(1)*q(3))];
        
        t = 1.1*RSEG*t;
        
        p = p + 1;
        
        isolated_seg_pos(:,p) = isolated_seg_pos(:,p-1) + t;
        
    end
    
    plot3((j-1)*L + isolated_seg_pos(1, 1 + (j-1)*NSEG: j*NSEG), isolated_seg_pos(2, 1 + (j-1)*NSEG: j*NSEG), isolated_seg_pos(3, 1 + (j-1)*NSEG: j*NSEG), 'b-');
 
    hold off;
    axis tight;
    axis equal;
    axis off;
    
end

sgtitle(sprintf('Scaled error = %g', min_error));
