close all
clear all
clc

addpath('/home/julien/ju_boulot/matlab/SHRU_programs/') %% local on moliere

sitename = 'SHRU1';
gps_site(1,:) = [72+54.4123/60 , -(159+1.0840/60)];
% depth(1,:)=[0 245 245.5 266.5]/100;
depth(1,:)=[0 2.5 5 7.5];


sitename = 'SHRU2';
gps_site(2,:) = [72+45.2347/60 , -(158+16.3243/60)];
% depth(2,:)=[0 260.5 266.5 250]/100;
depth(2,:)=[0 2.5 5 7.5];

sitename = 'SHRU3';
gps_site(3,:) = [72+40.6924/60 , -(157+54.6493/60)];
% depth(3,:)=[0 255.5, 280 252]/100;
depth(3,:)=[0 2.5 5 7.5];

sitename = 'SHRU4';
gps_site(4,:) = [72+36.6582/60 , -(157+32.2475/60)];
% depth(4,:)=[0 254.5, 274 247]/100;
depth(4,:)=[0 2.5 5 7.5];

sitename = 'SHRU5';
gps_site(5,:) = [72+54.4580/60 , -(157+29.2442/60)];
% depth(5,:)=[0 264, 518 782]/100;
depth(5,:)=[0 2.5 5 7.5];

figure
plot(gps_site(:,1), gps_site(:,2), 'o')

%% 2D
rr=2;
lat_ref=gps_site(rr,1);
lon_ref=gps_site(rr,2);

rot=0;
x_ref=0;
y_ref=0;
for rr=1:5
    lat=gps_site(rr,1);
    lon=gps_site(rr,2);
    [x(rr),y(rr)] = sub_transfer_LL_to_XY(lon,lat,lon_ref,lat_ref,rot,x_ref,y_ref);
end



figure
plot(x,y,'o')

pos_2d=[x; y; zeros(size(x))];
nor_2d=zeros(2,size(pos_2d,2)); 

%% 3D
pos_3d=[];
for rr=1:5
    lat=gps_site(rr,1);
    lon=gps_site(rr,2);
    [x(rr),y(rr)] = sub_transfer_LL_to_XY(lon,lat,lon_ref,lat_ref,rot,x_ref,y_ref);
    for zz=1:4    
        toto=[x(rr); y(rr); -depth(rr,zz)];
        pos_3d=[pos_3d toto];
    end
end
nor_3d=zeros(2,size(pos_3d,2)); 
figure
plot3(pos_3d(1,:),pos_3d(2,:),pos_3d(3,:),'.')

%%
% sensorArrayAnalyzer