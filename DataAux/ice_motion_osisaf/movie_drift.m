close all
clear all
clc

%%%%%% go to path /home/kfung/Downloads/CANAPE/DataAux/ice_motion_osisaf
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


ind=1;
for nn=1:size(nc_files,1)
    nameFile=nc_files(nn).name;
    pathFile=nc_files(nn).folder;
    
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


deltalon(deltalon>180)=deltalon(deltalon>180)-360;



%% Plot map with vector drift

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

clc

addpath('./quivermc_v4/');


%%% restrict to area of interest
toto=find(latitude<=latlimit(1) | latitude>=latlimit(2));
tata=find(longitude >= lonlimit_ok(1) & longitude <=lonlimit_ok(2));
tutu=union(toto, tata);



ind=1;
for nn=1:1:size(nc_files,1)
% for nn=100:10:130
    
    figure('visible','off');


    deltalat_ok=squeeze(deltalat(nn,:,:));
    deltalon_ok=squeeze(deltalon(nn,:,:));
    
    deltalat_ok(tutu)=NaN;
    deltalon_ok(tutu)=NaN;
    
       
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
    scale=2;
    quiverm(double(latitude),double(longitude),deltalat_ok,deltalon_ok,scale)
%     quivermc(double(latitude),double(longitude),deltalat_ok,deltalon_ok,scale, 'linewidth',2)
    title(datestr(date_num_ok(nn)))
    
    %%% add mooring
    plotm(gps_site(1), gps_site(2), 'xk','markersize',16,'linewidth',3)

    %%% save fig
%     h = get(0,'children');
%     scrsz = get(0,'ScreenSize');
%     set(h,'Position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)])  
%     set(gcf,'color','w'); 

    NameFig=['./fig_drift/' datestr(date_num_ok(nn),'yymmdd')];
    print(gcf,NameFig,'-dpng')
    
    close(gcf)
    disp([num2str(nn) '/' num2str(size(nc_files,1))])

end

disp('Done')