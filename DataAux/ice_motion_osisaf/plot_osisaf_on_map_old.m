close all
clear all
clc

plot_map=1;
write_csv_data=0;

sitename = 'SHRU1';
%%%% GPS data has longitude between -180 and 180
gps_site = [72.90687 , -159.01807];

nc_files = dir('./data/**/*.nc');

nn=105;
nameFile=nc_files(nn).name;
pathFile=nc_files(nn).folder;

finfo=ncinfo([pathFile '/' nameFile]);

latitude=ncread([pathFile '/' nameFile],'lat');
longitude=ncread([pathFile '/' nameFile],'lon');
%%%% Ice data has longitude between -180 and 180


dx=ncread([pathFile '/' nameFile],'dX');
dy=ncread([pathFile '/' nameFile],'dY');

d=sqrt(dx.^2+dy.^2);

ind_ok=[];
err=0.05;
err_step=0.025;
while isempty(ind_ok)    
    ind_lat = find(abs(gps_site(1)-latitude)<err);
    ind_lon = find(abs(gps_site(2)-longitude)<err);

    ind_ok=intersect(ind_lat,ind_lon);
    err=err+err_step;
end

if length(ind_ok)>1
    error('Put smaller err step in while loop')
end

[x,y]=ind2sub(size(latitude), ind_ok);


figure
pcolor(1:177,1:119,d)
shading flat
hold on
plot(y,x,'xk','markersize',16,'linewidth',3)


%% Better plot

figure
%%% set axes
latlimit=[65 90];
centralmeridian=-160;
maph=axesm('MapProjection','stereo','MapLatLimit',latlimit,  'Origin',[90 centralmeridian 0]);

%%% add grid
parallel=[70 75 80 85];
setm(maph,'Grid','on','Glinewidth',1,'PLineLocation',parallel)

%%% add grid labeling
MLabelLocation=[-120 -60 0 60 120 180];
Mpos=latlimit(1)-5;
PLabelMeridian=30;
setm(maph,'Fontangle','normal',...
  'FontSize',12,'fontweight','b',...
  'MeridianLabel','on',...
  'MLabelLocation',MLabelLocation,...
  'MLabelParallel',Mpos,...
  'ParallelLabel','on',...
  'PLabelLocation',parallel,...
  'PLabelMeridian',PLabelMeridian);

%%% add land
land = shaperead('landareas.shp', 'UseGeoCoords', true);
geoshow(maph, land, 'FaceColor',[0.80 0.80 0.80],'EdgeColor',0.30*[1 1 1]);

%%% plot data
surfm(double(latitude), double(longitude), d)

%%% add mooring
gps_site = [72.90687 , -159.01807];
hold on
plotm(gps_site(1),gps_site(2),'xk','markersize',16,'linewidth',3)

colorbar

%% make smaller map
min_lat=60;
min_lat_data=min(min(latitude));

ii=1;
Nx=119;
Ny=177;
while min_lat_data<min_lat
    x_small=ii:Nx+1-ii;
    y_small=ii:Ny+1-ii;
    latitude_small=latitude(x_small,y_small);
    min_lat_data=min(min(latitude_small));
    ii=ii+1;
end


% x_small=30:90;
% y_small=35:140;
% 
% latitude_small=latitude(x_small, y_small);
% longitude_small=longitude(x_small, y_small);
d_small=d(x_small, y_small);

figure
latlimit=[65 90];
centralmeridian=-160;
maph=axesm('MapProjection','stereo','MapLatLimit',latlimit,  'Origin',[90 centralmeridian 0]);

%%% add grid
parallel=[70 75 80 85];
setm(maph,'Grid','on','Glinewidth',1,'PLineLocation',parallel)

%%% add grid labeling
MLabelLocation=[-120 -60 0 60 120 180];
Mpos=latlimit(1)-5;
PLabelMeridian=30;
setm(maph,'Fontangle','normal',...
  'FontSize',12,'fontweight','b',...
  'MeridianLabel','on',...
  'MLabelLocation',MLabelLocation,...
  'MLabelParallel',Mpos,...
  'ParallelLabel','on',...
  'PLabelLocation',parallel,...
  'PLabelMeridian',PLabelMeridian);

%%% add land
land = shaperead('landareas.shp', 'UseGeoCoords', true);
geoshow(maph, land, 'FaceColor',[0.80 0.80 0.80],'EdgeColor',0.30*[1 1 1]);

%%% plot data
surfm(double(latitude), double(longitude), d)

%%% add mooring
gps_site = [72.90687 , -159.01807];
hold on
plotm(gps_site(1),gps_site(2),'xk','markersize',16,'linewidth',3)

colorbar
