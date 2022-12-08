clear

make_into_video = false;

wall = false;

FileName = 'ciliated_sphere_with_no_slip_surface';

if make_into_video
    
    video = VideoWriter(strcat(FileName,'_data_video.avi'));
    video.FrameRate = 10;
    
    open(video);
    
end

fid = fopen(strcat(FileName,'.dat'));
% fid_surf = fopen(strcat(FileName,'_surface.dat'));

% Read straight to the 2nd line of the filament data file to ignore comments
tline = fgetl(fid);
tline = fgetl(fid);

% tline_surf = fgetl(fid_surf);

% parameters = load(strcat(FileName,'.par'));
% 
% Nf = parameters(1);
% NperFil = parameters(2);
% Nsurf = parameters(3);
% 
% Nframes = floor(parameters(9)/10);
Nframes = 0;
Nf = 200;
NperFil = 20;

number_of_frames = 0;

while ischar(tline)
    
    number_of_frames = number_of_frames + 1;
        
        D = str2num(tline);
    
    for j=1:Nf
        
        FilD = D(2 + 13*NperFil*(j-1):1 + 13*NperFil*j);
        
        end_height = norm([FilD(1 + 13*(NperFil-1)),FilD(2 + 13*(NperFil-1)),FilD(3 + 13*(NperFil-1))]) - 200;
        
        scatter3(FilD(1),FilD(2),FilD(3),[],end_height/(2.2*NperFil),'filled');
        hold on;
        
    end
    
    caxis([0 1]);
    cbar = colorbar;
    cbar.Label.Interpreter = 'latex';
    cbar.Label.String = '$h/L$';
    
%     D = str2num(tline_surf);
%     
%     for j=1:Nsurf
%         
%         plot3(D(2+9*(j-1)),D(3+9*(j-1)),D(4+9*(j-1)),'k.');
%         
%     end
    
    %surf(200*s1,200*s2,200*s3,'FaceAlpha',1.0,'FaceColor',my_sphere_colour,'EdgeColor','none');
%     colormap winter;
%     surf(200*s1,200*s2,200*s3,'FaceAlpha',0.5);
%     shading interp;
    axis equal;
    axis([-250 250 -250 250 -250 250]);
    box on;
    grid off;
    
    hold off;
    
    if make_into_video
        
        writeVideo(video,getframe(gcf));
        
    else
        
        pause(0.05);
        
    end
    
    fprintf('Frame %i/%i drawn.\n',number_of_frames,Nframes);
    
    tline = fgetl(fid);
%     tline_surf = fgetl(fid_surf);
    
end

if make_into_video
    
    close(video);
    
end

fclose(fid);
% fclose(fid_surf);
