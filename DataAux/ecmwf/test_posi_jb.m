close all
clear all
clc

load 2017-06-01-00.mat

gps_site = [72.90687 , -159.01807];


longitude = longitude-180;

[c_lat,ind_lat] = min(abs(gps_site(1)-latitude));
[c_lon,ind_lon] = min(abs(gps_site(2)-longitude));

figure
imagesc(mwp)
hold on
plot(ind_lat,ind_lon,'xk','markersize',16,'linewidth',3)

clear longitude

load 2017-06-01-12.mat

longitude(longitude>180)=longitude(longitude>180)-360;


[c_lat,ind_lat] = min(abs(gps_site(1)-latitude));
[c_lon,ind_lon] = min(abs(gps_site(2)-longitude));
hold on
plot(ind_lat,ind_lon,'xr','markersize',16,'linewidth',3)