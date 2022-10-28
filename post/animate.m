function animate(sim_name,... % e.g. 'R_equals_1L_sphere6_f_220'
    start_fraction,... % number between 0 and 1 (inclusive)
    end_fraction,... % number between start_fraction and 1 (inclusive)
    show_figure,... % true or false
    view_vec,... % 2, 3 or a vector of the form [az el]
    show_time,... % true or false
    num_shown_frames,... % an integer >= 1
    body_fixed_view,... % true or false
    type,... % 'filament', 'tip_disp_vector_field' or 'prescribed'
    show_rotation_axis,... % true or false
    show_body_fixed_frame... % true or false
    )

% N.B. Only works for NBOD=1 because I don't currently run multi-body simulations.
%
% N.B. If the simulation was run in the normal git directory on a machine known as "remote" in your ssh
% config file, all strictly necessary data files can be downloaded by running the following bash command:
%
% for type in par reference state; do scp remote:~/phd_codes/tethered_filaments/sim_name*$type* .; done
%
% Some options, such as show_rotation_axis = true, may require additional files.

% -------------------------------------------------------------------- %
% Open the data files
% -------------------------------------------------------------------- %

fid_body = fopen([sim_name, '_body_states.dat']);
fid_seg = fopen([sim_name, '_seg_states.dat']);

par = read_parameter_file(sim_name);

if ~isfield(par, 'DL')
    
    par.DL = 2.2*RSEG; % This is the default value in old follower-force simulations.
    
end

if ~isfield(par, 'STEPS_PER_PERIOD')
    
    par.STEPS_PER_PERIOD = round(36.3833/DT); % Use the time scale for follower-force simulations. N.B. This time-scale is probably 'wrong' for any sim old enough to not save STEPS_PER_PERIOD anyway, because it almost certainly used pre-correction wall-RPY.
    
end

% This value is only setting axis limits etc. so it's fine that this is an
% overestimate for prescribed-motion filaments where we use L =
% DL*(NSEG-1). In fact, for beat patterns interpolated from real data the
% simulation L is only the average length over a period (in general) and
% having this margin of error is a good thing. Anywhere we want to use the
% actual physical scale from the simulations should call par.FIL_LENGTH.
L = par.DL * par.NSEG;

num_frames = 1 + floor((par.TOTAL_TIME_STEPS-1)/par.PLOT_FREQUENCY_IN_STEPS);

if exist([sim_name, '_blob_references.dat'], 'file') == 2

    is_planar_sim = false;

    if show_rotation_axis

        W = load([sim_name '_body_vels.dat']);
        W = W(:,5:7);

    end

else

    is_planar_sim = true;

    show_rotation_axis = false;

end

fil_references = reshape(load([sim_name, '_fil_references.dat']), 3, par.NFIL);

start_fraction = min(1, start_fraction);
start_fraction = max(0, start_fraction);
end_fraction = min(1, end_fraction);
end_fraction = max(end_fraction, start_fraction);

try

    % We may have not bothered to download the backup file, but if we
    % happen to have it then we can use it to get a tight upper bound on
    % how far the simulation reached.
    temp = load([sim_name '.backup']);
    goal_steps = par.TOTAL_TIME_STEPS;
    upper_bound_on_step_reached = min([temp(1) + par.PLOT_FREQUENCY_IN_STEPS, goal_steps]);

    if upper_bound_on_step_reached <= start_fraction*goal_steps

        % We expect to reach the end of the data before we actually start
        % plotting.
        fprintf('\n Backup file suggests the simulation ran for at most %g%% of the maximum steps.\n The choice start_fraction = %g is not expected to produce results. \n\n', 100*upper_bound_on_step_reached/par.TOTAL_TIME_STEPS, start_fraction);
        return;
        
    end
    
    % If we're asking our video to finish after the end of the file, reduce
    % the requested video length.
    end_fraction = min([end_fraction, upper_bound_on_step_reached/par.TOTAL_TIME_STEPS]);

end

% -------------------------------------------------------------------- %
% Prepare the video file
% -------------------------------------------------------------------- %

if body_fixed_view

    video_name = 'body_frame';

else

    video_name = 'lab_frame';

end

video_name = [sim_name '_' video_name '_' type '_animation_from_' num2str(start_fraction*par.TOTAL_TIME_STEPS/par.STEPS_PER_PERIOD) 'T_to_' num2str(end_fraction*par.TOTAL_TIME_STEPS/par.STEPS_PER_PERIOD) 'T'];

if show_rotation_axis

    video_name = [video_name '_with_rotation_axis'];

end

if show_body_fixed_frame

    video_name = [video_name '_with_body_fixed_frame'];

end

try

    % .mp4 is better than .avi in terms of space, but some remote machines
    % refuse to do .mp4
    video = VideoWriter(video_name, 'MPEG-4');

catch

    video = VideoWriter(video_name);

end

video.Quality = 100;
video.FrameRate = 30;
open(video);

% -------------------------------------------------------------------- %
% Initialise arrays etc.
% -------------------------------------------------------------------- %

read_line = @(fid, line_length) textscan(fid, repmat('%f', [1 line_length]), 1, 'CommentStyle', '%', 'Delimiter', ' ');

seg_pos = NaN(3, par.NSEG*par.NFIL + par.NFIL - 1, num_shown_frames);

if strcmp(type, 'tip_disp_vector_field')

    tip_disp = NaN(3, par.NFIL, num_shown_frames);
    base_pos = NaN(3, par.NFIL, num_shown_frames);

end

if is_planar_sim

    xmin = min(fil_references(1,:)) - L;
    xmax = max(fil_references(1,:)) + L;
    ymin = min(fil_references(2,:)) - L;
    ymax = max(fil_references(2,:)) + L;
    zmin = 0;
    zmax = 1.1*L;

else

    % TODO: How could I generalise this to non-spherical bodies?

    [s1, s2, s3] = sphere(100);

    sphere_radius = norm(fil_references(:,1));

    xmin = -sphere_radius - L;
    xmax = sphere_radius + L;
    ymin = -sphere_radius - L;
    ymax = sphere_radius + L;
    zmin = -sphere_radius - L;
    zmax = sphere_radius + L;

    sphere_radius = sphere_radius - 1.1*par.RSEG; % TODO: This is an old fixed difference and should be done based on reference files.
    s1 = sphere_radius*s1;
    s2 = sphere_radius*s2;
    s3 = sphere_radius*s3;

end

% -------------------------------------------------------------------- %
% Ensure the figure borders are minimal
% -------------------------------------------------------------------- %
% Note: the quality of this 'border removal' is dependent on the choice of
% view(...) and consequently how well the axis box fits the enclosed image.

% TO-DO: Is there a less 'hacky' way of doing this?

if show_figure

    h = figure;

else

    h = figure('visible', 'off');

end

h.CurrentAxes = axes;
ax = h.CurrentAxes;
axis(ax, 'equal');
axis(ax, [xmin xmax ymin ymax zmin zmax]);
view(ax, view_vec);
box(ax, 'on');
xticks(ax, []);
yticks(ax, []);
zticks(ax, []);
M = print(h, '-RGBImage', '-r1000');
M1 = double(M(:,:,1));
[row, col] = find(M1 - 255); % uint8 runs from 0 to 255
c = zeros(4,1);
r = zeros(4,1);
[c(1), id] = min(col); r(1) = row(id);
[c(2), id] = max(col); r(2) = row(id);
[r(3), id] = min(row); c(3) = col(id);
[r(4), id] = max(row); c(4) = col(id);
ax.Units = 'pixels';
h.Units = 'pixels';
height = max(r) - min(r);
width = max(c) - min(c);
if (height > width) % keep shortest dimension same as standard figure window
    ax.Position(4) = ax.Position(3)*height/width;
    h.Position(4) = h.Position(3)*height/width;
else
    ax.Position(3) = ax.Position(4)*width/height;
    h.Position(3) = h.Position(4)*width/height;
end
h.Position(1:2) = [100 100];
set(ax, 'Units', 'normalized', 'Position', [0 0 1 1]);
set(h, 'color', 'w'); % white background for better contrast

% -------------------------------------------------------------------- %
% Attempt to make the filament width accurate
% -------------------------------------------------------------------- %

[seg_s1, seg_s2, seg_s3] = sphere(100);
hold on;
seg_sphere_handle = surf(ax, par.RSEG*seg_s1 + 0.5*(xmin + xmax), par.RSEG*seg_s2 + 0.5*(ymin + ymax), par.RSEG*seg_s3 + 0.5*(zmin + zmax), 'EdgeColor', 'none', 'FaceColor', [0 0 0]);
axis off;
M = print(h, '-RGBImage', '-r1000');
M1 = double(M(:,:,1));
[row, col] = find(M1 - 255); % uint8 runs from 0 to 255
c = zeros(4,1);
r = zeros(4,1);
[c(1), id] = min(col); r(1) = row(id);
[c(2), id] = max(col); r(2) = row(id);
[r(3), id] = min(row); c(3) = col(id);
[r(4), id] = max(row); c(4) = col(id);
seg_height = max(r) - min(r);
seg_width = max(c) - min(c);
seg_size = 0.5*(seg_width + seg_height);
set(ax, 'Units', 'points');
line_width = 0.5*seg_size*(ax.Position(3)/size(M1,2) + ax.Position(4)/size(M1,1)); % take the average of the four possible widths we could calculate.
set(ax, 'Units', 'normalized', 'Position', [0 0 1 1]);
delete(seg_sphere_handle);

% -------------------------------------------------------------------- %
% Create the animation
% -------------------------------------------------------------------- %

hold on;

if ~is_planar_sim
    surf(ax, s1, s2, s3, 'Edgecolor', 'none', 'FaceColor', [0.8 0.8 0.8]);
end

if show_time
     h_string = text(ax, 0, 1, '$t/T = 0$', 'Interpreter', 'latex', 'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 24);
end

colours = gray(num_shown_frames+1);

for n = 1:num_frames

    try % fails if we reach the end of the file early; i.e. the simulation didn't run to the final time-step.

        D_body = read_line(fid_body, 8);
        
        if strcmp(type, 'prescribed')
            
            D_seg = read_line(fid_seg, 1 + par.NFIL*par.NSEG*3);
            
        else
            
            D_seg = read_line(fid_seg, 1 + par.NFIL*par.NSEG*4);
        
        end

        if ((n > start_fraction*num_frames) && (n <= end_fraction*num_frames))
            
            curr_time  = D_body{1}/par.STEPS_PER_PERIOD;
            
            D_body = cell2mat(D_body(2:end));

            if is_planar_sim

                R = 1;

            else

                q = D_body(4:7);
                qsq = q.^2;
                R = [1 - 2*(qsq(3) + qsq(4)), 2*(q(2)*q(3) - q(1)*q(4)), 2*(q(2)*q(4) + q(1)*q(3));...
                    2*(q(2)*q(3) + q(1)*q(4)), 1 - 2*(qsq(2) + qsq(4)), 2*(q(3)*q(4) - q(1)*q(2));...
                    2*(q(2)*q(4) - q(1)*q(3)), 2*(q(3)*q(4) + q(1)*q(2)), 1 - 2*(qsq(2) + qsq(3))];

            end

            D_seg = cell2mat(D_seg(2:end));

            try

                delete(handles);

                if show_rotation_axis

                    delete(rot_axis_handle);

                end
                
                if show_body_fixed_frame
                    
                    delete(frame_handle);
                    
                end

            end

            for m=1:num_shown_frames-1

                seg_pos(:,:,m) = seg_pos(:,:,m+1);

                if strcmp(type, 'tip_disp_vector_field')

                    tip_disp(:,:,m) = tip_disp(:,:,m+1);
                    base_pos(:,:,m) = base_pos(:,:,m+1);

                end

            end

            for j=par.NFIL:-1:1
                
                if strcmp(type, 'prescribed')
                    
                    Dfil = D_seg(1 + 3*(j-1)*par.NSEG: 3*j*par.NSEG);
                    
                else

                    Dfil = D_seg(1 + 4*(j-1)*par.NSEG: 4*j*par.NSEG);
                
                end

                p = j + (j-1)*par.NSEG;

                if body_fixed_view

                    seg_pos(:,p,end) = fil_references(:,j);

                else

                    seg_pos(:,p,end) = R*fil_references(:,j);

                end

                if strcmp(type, 'tip_disp_vector_field')

                    base_pos(:,j,end) = seg_pos(:,p,end);

                end

                for k=2:par.NSEG
                    
                    if strcmp(type, 'prescribed')
                        
                        id = 3*(k-1) + 1;
                        p = p + 1;
                        
                        if body_fixed_view
                            
                            seg_pos(:,p,end) = R'*(Dfil(id:id+2)' - D_body(1:3)');
                            
                        else
                            
                            seg_pos(:,p,end) = Dfil(id:id+2)' - D_body(1:3)';
                            
                        end
                        
                    else

                        id = 4*(k-2) + 1;

                        q = Dfil(id:id+3);
                        qsq = q.^2;
                        t = [1 - 2*(qsq(3) + qsq(4)); 2*(q(2)*q(3) + q(1)*q(4)); 2*(q(2)*q(4) - q(1)*q(3))];

                        id = id + 4;

                        q = Dfil(id:id+3);
                        qsq = q.^2;
                        t = t + [1 - 2*(qsq(3) + qsq(4)); 2*(q(2)*q(3) + q(1)*q(4)); 2*(q(2)*q(4) - q(1)*q(3))];

                        t = 0.5*par.DL*t;

                        p = p + 1;

                        if body_fixed_view

                            seg_pos(:,p,end) = seg_pos(:,p-1,end) + R'*t;

                        else

                            seg_pos(:,p,end) = seg_pos(:,p-1,end) + t;

                        end
                    
                    end

                end

                if strcmp(type, 'tip_disp_vector_field')

                    tip_disp(:,j,end) = seg_pos(:,p,end) - base_pos(:,j,end);

                    if ~is_planar_sim % project into tangent plane if on a spherical surface

                        tip_disp(:,j,end) = tip_disp(:,j,end) - dot(tip_disp(:,j,end), base_pos(:,j,end))*base_pos(:,j,end)/dot(base_pos(:,j,end), base_pos(:,j,end));

                    end

                end

            end

            for m=num_shown_frames:-1:1

                if (strcmp(type, 'filament') || strcmp(type, 'prescribed'))

                    handles(m) = plot3(ax, seg_pos(1,:,m), seg_pos(2,:,m), seg_pos(3,:,m), '-', 'LineWidth', line_width, 'Color', colours(m,:));

                elseif strcmp(type, 'tip_disp_vector_field')

                    if is_planar_sim

                        handles(m) = quiver(ax, base_pos(1,:,m), base_pos(2,:,m), tip_disp(1,:,m), tip_disp(2,:,m), 'Color', colours(m,:));

                    else

                        handles(m) = quiver3(ax, base_pos(1,:,m), base_pos(2,:,m), base_pos(3,:,m), tip_disp(1,:,m), tip_disp(2,:,m), tip_disp(3,:,m), 'Color', colours(m,:));

                    end

                end

            end

            if show_rotation_axis

                Wvec = W(n,:)';
                Wvec = (L + sphere_radius)*Wvec/norm(Wvec);

                if body_fixed_view

                    Wvec = R' * Wvec;

                end

                rot_axis_handle = quiver3(ax, [0 0], [0 0], [0 0], Wvec(1)*[1 -1], Wvec(2)*[1 -1], Wvec(3)*[1 -1], 0, 'r-');

            end
            
            if show_body_fixed_frame

                if body_fixed_view

                    frame_handle = quiver3(ax, [0 0 0], [0 0 0], [0 0 0], (L + sphere_radius)*[1 0 0], (L + sphere_radius)*[0 1 0], (L + sphere_radius)*[0 0 1], 0, 'b-');
                    
                else
                    
                    frame_handle = quiver3(ax, [0 0 0], [0 0 0], [0 0 0], (L + sphere_radius)*R(1,:), (L + sphere_radius)*R(2,:), (L + sphere_radius)*R(3,:), 0, 'b-');

                end

            end

            if show_time
                h_string.String = sprintf('$t/T = %g$', curr_time);
            end

            if show_figure

                frame = getframe(h);

            else

                % Bypass the need to use getframe(), which sometimes fails
                % when 'visible' is set to 'off'.
                print(h, 'temp_frame_for_animation', '-djpeg', '-r300');
                frame = imread('temp_frame_for_animation.jpg');

            end

            writeVideo(video, frame);
            fprintf('Written frame %i.\n', n);
            
        elseif (n > end_fraction*num_frames)
            
            break;

        end

    catch my_err
        
        disp(my_err)

        close(video);

        if ~show_figure
            close(h);
        end

        fclose all;
        return;

    end

end

close(video);

if ~show_figure
    close(h);
end

fclose all;

end
