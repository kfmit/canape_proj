close all
clear all
clc


nc_files = dir('./data/**/*.nc');

% sitename = 'SHRU1';
% gps_site = [72+54.4123/60 , -(159+1.0840/60)];

% sitename = 'SHRU2';
% gps_site = [72+45.2347/60 , -(158+16.3243/60)];
% 
% sitename = 'SHRU3';
% gps_site = [72+40.6924/60 , -(157+54.6493/60)];
% 
% sitename = 'SHRU4';
% gps_site = [72+36.6582/60 , -(157+32.2475/60)];
% 
sitename = 'SHRU5';
gps_site = [72+54.4580/60 , -(157+29.2442/60)];


nn=105;
nameFile=nc_files(nn).name;
pathFile=nc_files(nn).folder;

latitude=ncread([pathFile '/' nameFile],'lat');
longitude=ncread([pathFile '/' nameFile],'lon');
%%%% Ice data has longitude between -180 and 180

%% Restrict to one day a month
n_ok=[];
for nn=1:size(nc_files,1)
    nameFile=nc_files(nn).name;
    if str2num(nameFile(end-8:end-7)) == 2
        n_ok=[n_ok nn];
%         nameFile
    end
end

deb=3; 
ind=1;
for nn=deb:deb+5
    nameFile=nc_files(n_ok(nn)).name;
    pathFile=nc_files(n_ok(nn)).folder;
    
    date_vec=nameFile(end-14:end-3);
    date_num=datenum(date_vec,'yyyymmddHHMM');
    date_num_ok(ind)=date_num-1;
    date_vec_ok=datestr(date_num_ok(ind),'yyyymmddHHMMSS');
    datevect{ind,1} = date_vec_ok;
    
    dx=ncread([pathFile '/' nameFile],'dX');
    dy=ncread([pathFile '/' nameFile],'dY');
%     lat0=ncread([pathFile '/' nameFile],'lat0');
%     lon0=ncread([pathFile '/' nameFile],'lon0');
    lat1=ncread([pathFile '/' nameFile],'lat1');
    lon1=ncread([pathFile '/' nameFile],'lon1');
    d(ind,:,:)=sqrt(dx.^2+dy.^2);
    deltalat(ind,:,:)=lat1-latitude;
    deltalon(ind,:,:)=lon1-longitude;
    
    ind=ind+1;
end

% %%%% May 15 is only Nan
% nn=deb+5;
% nameFile=nc_files(n_ok(nn)-14).name;
% pathFile=nc_files(n_ok(nn)-14).folder;
% 
% date_vec=nameFile(end-14:end-3);
% date_num=datenum(date_vec,'yyyymmddHHMM');
% date_num_ok(ind)=date_num-1;
% date_vec_ok=datestr(date_num_ok(ind),'yyyymmddHHMMSS');
% datevect{ind,1} = date_vec_ok;
% 
% dx=ncread([pathFile '/' nameFile],'dX');
% dy=ncread([pathFile '/' nameFile],'dY');
% %     lat0=ncread([pathFile '/' nameFile],'lat0');
% %     lon0=ncread([pathFile '/' nameFile],'lon0');
% lat1=ncread([pathFile '/' nameFile],'lat1');
% lon1=ncread([pathFile '/' nameFile],'lon1');
% d(ind,:,:)=sqrt(dx.^2+dy.^2);
% deltalat(ind,:,:)=lat1-latitude;
% deltalon(ind,:,:)=lon1-longitude;


deltalon(deltalon>180)=deltalon(deltalon>180)-360;
%% Plot map with drift amplitude


latlimit=[65 90];
centralmeridian=-160;
parallel=[70 75 80 85];
MLabelLocation=[-120 -60 0 60 120 180];
Mpos=latlimit(1)-5;
PLabelMeridian=30;
land = shaperead('landareas.shp', 'UseGeoCoords', true);

figure
for nn=1:6
    subplot(2,3,nn)

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
    surfm(double(latitude),double(longitude),squeeze(d(nn,:,:)))
    title(datestr(date_num_ok(nn)))
    
    c = colorbar;
    c.Label.String = 'Ice drift (km)';
    caxis([0 70])
    
    %%% add mooring
    plotm(gps_site(1), gps_site(2), 'xk','markersize',16,'linewidth',3)

end

%% Plot map with vector drift

latitude_limit_plot=80;
toto=find(latitude>latitude_limit_plot);

figure
for nn=1:6
    subplot(2,3,nn)

    deltalat_ok=squeeze(deltalat(nn,:,:));
    deltalon_ok=squeeze(deltalon(nn,:,:));
    
    deltalat_ok(toto)=NaN;
    deltalon_ok(toto)=NaN;
    
    
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
    quiverm(double(latitude),double(longitude),deltalat_ok,deltalon_ok,scale)
    title(datestr(date_num_ok(nn)))
    
    %%% add mooring
    plotm(gps_site(1), gps_site(2), 'xk','markersize',16,'linewidth',3)
end