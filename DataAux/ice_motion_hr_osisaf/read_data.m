close all
clear all
clc

nc_files = dir('./data_nc/*.nc');

sitename = 'SHRU5';
gps_site = [72+54.4580/60 , -(157+29.2442/60)];

nn=450;
nameFile=nc_files(nn).name;
pathFile=nc_files(nn).folder;

latitude=ncread([pathFile '/' nameFile],'lat');
longitude=ncread([pathFile '/' nameFile],'lon');

toto=ncinfo([pathFile '/' nameFile]);
ttt=toto.Attributes(22).Value

dx=ncread([pathFile '/' nameFile],'dX');
dy=ncread([pathFile '/' nameFile],'dY');
dx(dx<-900)=NaN;
dy(dy<-900)=NaN;
d=sqrt(dx.^2+dy.^2);

lat1=ncread([pathFile '/' nameFile],'lat1');
lon1=ncread([pathFile '/' nameFile],'lon1');
deltalat(:,:)=lat1-latitude;
deltalon(:,:)=lon1-longitude;
deltalat(deltalat<-900)=NaN;
deltalon(deltalon<-900)=NaN;
deltalon(abs(deltalon)>100)=NaN;

latlimit=[65 85];
dlon=40;
lonlimit=[gps_site(2)-dlon gps_site(2)+dlon];
lonlimit_ok=[lonlimit(2) lonlimit(1)+360];

centralmeridian=-160;
parallel=[70 75 80];
MLabelLocation=[-180:20:180];
Mpos=latlimit(1)+2;
PLabelMeridian=30;
land = shaperead('landareas.shp', 'UseGeoCoords', true);

figure
%%% set axes
maph=axesm('MapProjection','lambertstd','MapLatLimit',latlimit,'MapLonLimit',lonlimit);

%%% add grid
setm(maph,'Grid','on','Glinewidth',1,'PLineLocation',parallel, 'MLineLocation',MLabelLocation)

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
surfm(double(latitude),double(longitude),squeeze(d(:,:)))

c = colorbar;
c.Label.String = 'Ice drift (km)';
% caxis([0 70])

%%% add mooring
plotm(gps_site(1), gps_site(2), 'xk','markersize',16,'linewidth',3)

%%
figure
%%% set axes
maph=axesm('MapProjection','lambertstd','MapLatLimit',latlimit,'MapLonLimit',lonlimit);

%%% add grid
setm(maph,'Grid','on','Glinewidth',1,'PLineLocation',parallel, 'MLineLocation',MLabelLocation)

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
scale=200;
quiverm(double(latitude),double(longitude),deltalat,deltalon,scale)

c = colorbar;
c.Label.String = 'Ice drift (km)';
% caxis([0 70])

%%% add mooring
plotm(gps_site(1), gps_site(2), 'xk','markersize',16,'linewidth',3)    


%%
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
scale=2;
quiverm(double(latitude),double(longitude),deltalat,deltalon,scale)

c = colorbar;
c.Label.String = 'Ice drift (km)';
% caxis([0 70])

%%% add mooring
plotm(gps_site(1), gps_site(2), 'xk','markersize',16,'linewidth',3)   

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
% geoshow(maph, land, 'FaceColor',[0.80 0.80 0.80],'EdgeColor',0.30*[1 1 1]);
% 
% %%% add ice data
% surfm(double(latitude),double(longitude),edges)
% % title(datestr(datevect_num_ok(nn)))
% 
% %%% add mooring
% plotm(gps_site(1), gps_site(2), 'xk','markersize',16,'linewidth',3)
   