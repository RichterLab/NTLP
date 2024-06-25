clear
clc
close all

%The "history file" contains quantities that have been averaged either volumetrically or horizontally
%The frequency that this is written out is governed by the "ihist" paramater in params.in
%To see the contents of the file within Matlab, use the "ncdisp" command

%"history.nc" is a netCDF file which can be read using standard packages


fname = './history.nc';

%% Plot the volume-averaged relative humidity versus time

%Load the time and meanRH variables:
time = ncread(fname,'time');
meanRH = ncread(fname,'meanRH');

figure(1)
plot(time,meanRH,'linewidth',2)
xlabel('time [s]')
ylabel('Volume-average relative humidity [%]')


%% Plot the average radius of all particles versus time

radavg = ncread(fname,'radavg');

figure(2)
plot(time,radavg,'linewidth',2)
xlabel('time [s]')
ylabel('Average radius of all particles [m]')


%% Plot the water vapor mixing ratio versus time

txym = ncread(fname,'txym');  %"txym" contains BOTH potential temperature (theta) and water vapor mixing ratio (qv)
qxym = squeeze(txym(:,2,:));

zu = ncread(fname,'zu'); %"zu" is the vertical grid points at the locations which contain u, v, and scalars. "zw" is the height which contains w and fluxes


%Multiple ways of plotting the vertical profile versus height

%1. Animate

Nt = length(time); %Number of time steps
stride = 10; %Don't necessarily want to animate each frame
figure(3)
for i=1:stride:Nt

    figure(3); clf
    plot(qxym(:,i),zu(:),'linewidth',2)
    xlabel('qv [kg/kg]')
    ylabel('z [m]')
    drawnow

end


%%
%2. Plot as a time-height contour
Nz = length(zu(:,1)); %Number of u-level grid points in z

%Make 2D arrays to plot agains
Z = repmat(zu(:,1),1,Nt);
T = repmat(time(:)',Nz,1);

clevels = linspace(0.01,0.02,100);

figure(4)
contourf(T,Z,qxym,clevels,'edgecolor','none')
colorbar
title('qv [kg/kg]')
xlabel('time [s]')
ylabel('z [m]')


%3. Perform a time average

%Need to know the time index beyond which it is appropriate to average (after statistical stationarity)
%Here take 100

tidx = 100;

qxym_avg = mean(qxym(:,tidx:end),2);

figure(5)
plot(qxym_avg,zu(:,1),'linewidth',2)
xlabel('<qv> [kg/kg]')
ylabel('z [m]')
