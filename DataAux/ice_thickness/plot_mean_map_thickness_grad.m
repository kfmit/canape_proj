close all
clear all
clc

nc_files = dir('./smos_sea_ice_thickness/*.nc');

% sitename = 'SHRU1';
% gps_site = [72+54.4123/60 , -(159+1.0840/60)];
% 
% sitename = 'SHRU2';
% gps_site = [72+45.2347/60 , -(158+16.3243/60)];

% sitename = 'SHRU3';
% gps_site = [72+40.6924/60 , -(157+54.6493/60)];
% 
% sitename = 'SHRU4';
% gps_site = [72+36.6582/60 , -(157+32.2475/60)];
% 
sitename = 'SHRU5';
gps_site = [72+54.4580/60 , -(157+29.2442/60)];

latitude = ncread(['./smos_sea_ice_thickness/' nc_files(1).name], 'latitude');
longitude = ncread(['./smos_sea_ice_thickness/' nc_files(1).name], 'longitude');


%% Restrict to one day a month
n_ok=[];
for nn=1:size(nc_files,1)
    nameFile=nc_files(nn).name;
    if str2num(nameFile(end-4:end-3)) == 1
        n_ok=[n_ok nn];
    end
end


%%

ave_vec=-15:14;
for nn=2:length(n_ok)

    thick_temp=zeros(size(latitude));
    T_temp=zeros(size(latitude));
    for mm=ave_vec
        nameFile=nc_files(n_ok(nn)+mm).name;
        datevect{nn,1} = [nameFile(end-10:end-3),'000000'];   
        datevect_num(nn)=datenum([str2num(datevect{nn,1}(1:4)),str2num(datevect{nn,1}(5:6)),str2num(datevect{nn,1}(7:8)),0,0,0]);

        thick_temp=thick_temp+ncread(['./smos_sea_ice_thickness/' nameFile], 'sea_ice_thickness');
        T_temp=T_temp+ncread(['./smos_sea_ice_thickness/' nameFile], 'Tsurf');

    end

    sea_ice_thickness{nn} = thick_temp / length(ave_vec);
    Tsurf{nn} = T_temp / length(ave_vec);

end


%% Plot Thickness



latlimit=[65 80];
dlon=40;
lonlimit=[gps_site(2)-dlon gps_site(2)+dlon];
centralmeridian=-160;
parallel=[70 75 80];
MLabelLocation=[-180:20:180];
Mpos=latlimit(1)+2;
PLabelMeridian=30;
land = shaperead('landareas.shp', 'UseGeoCoords', true);


figure

for nn=3:7
    subplot(2,3,nn-2)
    thickness=sea_ice_thickness{nn};
    
    thickness(thickness<0)=NaN;
    [Gmag,Gdir] = imgradient(thickness);

    %%% set axes
%     maph=axesm('MapProjection','stereo','MapLatLimit',latlimit,'MapLonLimit',lonlimit, 'Origin',[90 centralmeridian 0]);
%     maph=axesm('MapProjection','mercator','MapLatLimit',latlimit,'MapLonLimit',lonlimit);
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
    surfm(latitude,longitude,Gmag)
    title([datestr(datevect_num(nn)) ' \pm 15 days'] )
    
    c = colorbar;
    c.Label.String = 'Ice Thickness gradient';
    caxis([0 1.5])

    
    
        
    %%% add mooring
    plotm(gps_site(1), gps_site(2), 'xk','markersize',16,'linewidth',3)
end

