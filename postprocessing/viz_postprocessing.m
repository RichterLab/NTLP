tic

clear
clc
close all
%%

%The "viz file" contains slices of velocities and scalars
%The frequency that this is written out is governed by the "i_viz" paramater in params.in
%To see the contents of the file within Matlab, use the "ncdisp" command
%By default the slices are at the midplanes in each direction; to change, see the subroutine "write_viz_netcdf"

%"viz.nc" is a netCDF file which can be read using standard packages

% I modified this script to save each of the animations at 32 frames per second as an mp4 file in the folder that the .m file is run in. To change
% the framerate of all animations (and therefore the length as the number of frames is constant) change the variable called "fps" just below. To
% change the framerate of an individual animation, change the value of vidfile.FrameRate after the animation you want to affect to the desired
% framerate. DO NOT change the size of the figure while the animation is running as it will not write the .mp4 file. You can change the resolution
% in the figure object to change the resolution of the animation (it will be very large by default and may not show up on screen if your monitor is
% small)
% 
% - Ben Roper 2024

fps = 32;

fname = './viz.nc';

cmap = colormap("turbo"); % This makes the matlab native plots have more contrast


%% Animate an xy-slice of potential temperature
tic
time = ncread(fname,'time');
Nt = length(time);

x = ncread(fname,'x');
Nx = length(x);

y = ncread(fname,'y');
Ny = length(y);


X = repmat(x,1,Ny);
Y = repmat(y',Nx,1);

t_xy = ncread(fname,'t_xy');

min_t = min(min(min(t_xy)));
max_t = max(max(max(t_xy)));
clevels = linspace(min_t,max_t,100);

fig = figure('units','pixels','position',[0 0 1440 1080],'Colormap',cmap);
for i=1:Nt
    clf
    contourf(X,Y,t_xy(:,:,i),clevels,'edgecolor','none')
    colorbar
    clim([min_t max_t])
    xlabel('x [m]')
    ylabel('y [m]')
    title('XY Slice of Potential Temperature')
    axis equal
    movie1(i) = getframe(fig);
end

vidfile = VideoWriter('XY Slice of Potential Temperature 32fps','MPEG-4');
vidfile.FrameRate = fps;

open(vidfile)
writeVideo(vidfile,movie1)
close(vidfile)
close all
toc
%% Animate yz-slice of vertical velocity [NOTE: truncates the z = 0 point to make the same size as the u, v, scalar arrays]
tic
z = ncread(fname,'zu');
Nz = length(z);

w_yz = ncread(fname,'w_yz');

Y = repmat(y,1,Nz);
Z = repmat(z',Ny,1);

min_w = min(min(min(w_yz)));
max_w = max(max(max(w_yz)));
clevels = linspace(min_w,max_w,100);

fig2 = figure('units','pixels','position',[0 0 1440 1080],'Colormap',cmap);

for i=1:Nt
    clf
    contourf(Y,Z,w_yz(:,:,i),clevels,'edgecolor','none')
    colorbar
    clim([min_w max_w])
    xlabel('y [m]')
    ylabel('z [m]')
    title('YZ Slice of Vertical Velocity')
    axis equal
    movie2(i) = getframe(fig2);
end

vidfile = VideoWriter('YZ Slice of Vertical Velocity 32fps','MPEG-4');
vidfile.FrameRate = fps;

open(vidfile)
writeVideo(vidfile,movie2)
close(vidfile)

close all
toc

%% QUICK RENDERS

% the following two animations will be pixelated and square, but will
% animate much faster. To use this part of the script, the add on titled
% "SC - powerful image rendering" must be installed. After install,
% uncomment the code below and run the sections
%
% - Ben Roper 2024

% %% Quick render of w_yz
% tic
% 
% figure
% for ii = 1:Nt
%     sc(  flipud((w_yz(:,:,ii))')  ,[min_w max_w],'cold');
%     movie3(ii) = getframe;
% end
% 
% vidfile = VideoWriter('YZ_W 24 fps','MPEG-4');
% vidfile.FrameRate = fps;
% 
% open(vidfile)
% writeVideo(vidfile,movie3)
% close(vidfile)
% close all
% 
% toc
% %% Quick render of t_xy
% tic
% 
% figure
% for ii = 1:Nt
%     sc(  flipud((t_xy(:,:,ii))')  ,[min_t max_t],'cold');
%     movie4(ii) = getframe;
% end
% 
% vidfile = VideoWriter('XY_T 24 fps','MPEG-4');
% vidfile.FrameRate = fps;
% 
% open(vidfile)
% writeVideo(vidfile,movie4)
% close(vidfile)
% close all
% 
% toc

%%
toc