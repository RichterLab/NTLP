clear; clc; close all;

% This script is just the quick visualizations from viz_postprocessing.m
% including a separate animation of the velocity magnitude on the
% same yz plane as the vertical velocity.
% This script requires the addon "SC - powerful image rendering"
%
% - Ben Roper 2024

fps = 32;
fname = './viz.nc';
time = ncread(fname,'time');
Nt = length(time);
%% T_xy
t_xy = ncread(fname,'t_xy');
min_t = min(min(min(t_xy)));
max_t = max(max(max(t_xy)));
for ii = 1:Nt
    sc(flipud((t_xy(:,:,ii))'),[min_t max_t],'cold');
    movie1(ii) = getframe;
end
vidfile = VideoWriter('XY_T quick','MPEG-4');
vidfile.FrameRate = fps;
open(vidfile)
writeVideo(vidfile,movie1)
close(vidfile)
close all
%% W_yz
w_yz = ncread(fname,'w_yz');
min_w = min(min(min(w_yz)));
max_w = max(max(max(w_yz)));
for ii = 1:Nt
    sc(flipud((w_yz(:,:,ii))'),[min_w max_w],'cold');
    movie2(ii) = getframe;
end
vidfile = VideoWriter('YZ_W quick','MPEG-4');
vidfile.FrameRate = fps;
open(vidfile)
writeVideo(vidfile,movie2)
close(vidfile)
close all
%% vmag_yz
u_yz = ncread(fname,'u_yz');
v_yz = ncread(fname,'v_yz');
w_yz = ncread(fname,'w_yz');
vmag = zeros(size(u_yz));
tic
for ii = 1:Nt
    vmag(:,:,ii) = sqrt((u_yz(:,:,ii).^2) + (v_yz(:,:,ii).^2) + (w_yz(:,:,ii).^2));
    min_vmag = min(min(min(vmag)));
    max_vmag = max(max(max(vmag)));
    sc((vmag(:,:,ii)'),[min_vmag max_vmag],'hot')
    movie3(ii) = getframe; 
end
toc
vidfile = VideoWriter('vmag quick 32 fps','MPEG-4');
vidfile.FrameRate = fps;
open(vidfile)
writeVideo(vidfile,movie3)
close(vidfile)
close all