close all
clear all
clc



nc_files = dir('./data_nc/*.nc');

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


nn=1;
nameFile=nc_files(nn).name;

latitude = ncread(['./data_nc/' nc_files(1).name], 'latitude');
%%%% Satellite map has longitude between 360 and 0
longitude = ncread(['./data_nc/' nc_files(1).name], 'longitude');
% longitude(longitude>180)=longitude(longitude>180)-360;



%% Restrict to one day a month
n_ok=[];
for nn=1:size(nc_files,1)
    nameFile=nc_files(nn).name;
    if str2num(nameFile(7:8)) == 1
        n_ok=[n_ok nn];
    end
end
    
ave_vec=-15:15;
for nn=1:6
    ice_frac_ok{nn}=zeros(size(latitude));
    for mm=ave_vec
        nameFile=nc_files(n_ok(nn+2)+mm).name;
        ice_frac_ok{nn} = ice_frac_ok{nn} + ncread(['./data_nc/' nameFile], 'sea_ice_area_fraction');
    end
    ice_frac_ok{nn}=ice_frac_ok{nn}/length(ave_vec);

    nameFile=nc_files(n_ok(nn+2)).name;
    datevect{nn,1} = [nameFile(1:8),'000000'];
    datevect_num_ok(nn)=datenum([str2num(datevect{nn,1}(1:4)),str2num(datevect{nn,1}(5:6)),str2num(datevect{nn,1}(7:8)),0,0,0]);  
    
end

%% Plot map


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
    frac=ice_frac_ok{nn};
    
	frac(frac>100)=NaN;

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
    surfm(double(latitude),double(longitude),frac)
    title([datestr(datevect_num_ok(nn)) ' \pm 15 days'] )
     
    c = colorbar;
    c.Label.String = 'Ice fraction (%)';
    caxis([80 100])

    %%% add mooring
    plotm(gps_site(1), gps_site(2), 'xk','markersize',16,'linewidth',3)
end

figure
for nn=1:6
    subplot(2,3,nn)
    frac=ice_frac_ok{nn};
    
	frac(frac>100)=NaN;

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
    surfm(double(latitude),double(longitude),frac)
    title([datestr(datevect_num_ok(nn)) ' \pm 15 days'] )
     
    c = colorbar;
    c.Label.String = 'Ice fraction (%)';
    caxis([0 100])

    %%% add mooring
    plotm(gps_site(1), gps_site(2), 'xk','markersize',16,'linewidth',3)
end
