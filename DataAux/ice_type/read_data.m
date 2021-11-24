close all
clear all
clc



sitename = 'SHRU5';
gps_site = [72+54.4580/60 , -(157+29.2442/60)];

nc_files = dir('./data_nc/*.nc');
% nc_files = dir('./*.nc');

nn=31;
nameFile=[nc_files(nn).folder,filesep, nc_files(nn).name];

% finfo=ncinfo(nameFile)

latitude = ncread(nameFile, 'lat');
longitude = ncread(nameFile, 'lon');

ice_type = ncread(nameFile, 'ice_type');

types=imgradient(ice_type);
types(types>.5)=10;


toto=find(types==10);
latitude_types=double(latitude(toto));
longitude_types=double(longitude(toto));

%% Plot map

latlimit=[65 90];
centralmeridian=-160;
parallel=[70 75 80 85];
MLabelLocation=[-120 -60 0 60 120 180];
Mpos=latlimit(1)-5;
PLabelMeridian=30;
land = shaperead('landareas.shp', 'UseGeoCoords', true);

figure
%%% set axes
maph=axesm('MapProjection','stereo','MapLatLimit',latlimit,  'Origin',[90 centralmeridian 0]);

%%% add grid
setm(maph,'Grid','on','Glinewidth',1,'PLineLocation',parallel)

%%% add grid labeling
setm(maph,'Fontangle','normal',...
  'FontSize',12,'fontweight','b',...
  'MeridianLabel','on',...
  'MLabelLocation',MLabelLocation,...
  'MLabelParallel',Mpos,...
  'ParallelLabel','on',...
  'PLabelLocation',parallel,...
  'PLabelMeridian',PLabelMeridian);

%%% add land
geoshow(maph, land, 'FaceColor',[0.80 0.80 0.80],'EdgeColor',0.30*[1 1 1]);

%%% add ice data
surfm(double(latitude),double(longitude),ice_type)
% title(datestr(datevect_num_ok(nn)))

%%% add mooring
plotm(gps_site(1), gps_site(2), 'xk','markersize',16,'linewidth',3)

%%% add types
plotm(latitude_types, longitude_types, '.k','markersize',5,'linewidth',1)

% figure
% %%% set axes
% maph=axesm('MapProjection','stereo','MapLatLimit',latlimit,  'Origin',[90 centralmeridian 0]);
% 
% %%% add grid
% setm(maph,'Grid','on','Glinewidth',1,'PLineLocation',parallel)
% 
% %%% add grid labeling
% setm(maph,'Fontangle','normal',...
%   'FontSize',12,'fontweight','b',...
%   'MeridianLabel','on',...
%   'MLabelLocation',MLabelLocation,...
%   'MLabelParallel',Mpos,...
%   'ParallelLabel','on',...
%   'PLabelLocation',parallel,...
%   'PLabelMeridian',PLabelMeridian);
% 
% %%% add land
% geoshow(maph, land, 'FaceColor',[0.80 0.80 0.80],'typeColor',0.30*[1 1 1]);
% 
% %%% add ice data
% surfm(double(latitude),double(longitude),types)
% % title(datestr(datevect_num_ok(nn)))
% 
% %%% add mooring
% plotm(gps_site(1), gps_site(2), 'xk','markersize',16,'linewidth',3)