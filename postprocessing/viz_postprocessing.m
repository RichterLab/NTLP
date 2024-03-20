

%The "viz file" contains slices of velocities and scalars
%The frequency that this is written out is governed by the "i_viz" paramater in params.in
%To see the contents of the file within Matlab, use the "ncdisp" command
%By default the slices are at the midplanes in each direction; to change, see the subroutine "write_viz_netcdf"

%"viz.nc" is a netCDF file which can be read using standard packages


clear
clc
close all

fname = './viz.nc';

%% Animate an xy-slice of potential temperature

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

figure(1)
for i=1:Nt

    figure(1); clf
    contourf(X,Y,t_xy(:,:,i),clevels,'edgecolor','none')
    colorbar
    caxis([min_t max_t])
    xlabel('x [m]')
    ylabel('y [m]')
    drawnow

end


%% Animate yz-slice of vertical velocity [NOTE: truncates the z = 0 point to make the same size as the u, v, scalar arrays]

z = ncread(fname,'zu');
Nz = length(z);

w_yz = ncread(fname,'w_yz');

Y = repmat(y,1,Nz);
Z = repmat(z',Ny,1);

min_w = min(min(min(w_yz)));
max_w = max(max(max(w_yz)));
clevels = linspace(min_w,max_w,100);

figure(2)
for i=1:Nt

    figure(2); clf
    contourf(Y,Z,w_yz(:,:,i),clevels,'edgecolor','none')
    colorbar
    caxis([min_w max_w])
    xlabel('y [m]')
    ylabel('z [m]')
    drawnow

end
