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



nn=1;
nameFile=nc_files(nn).name;
land=ncread(['./smos_sea_ice_thickness/' nameFile], 'land');
ind_land=find(land);

imAlpha=ones(size(land));
imAlpha(ind_land)=0;



%%



datevect={};
vecice_frac=zeros(size(nc_files,1),4);

for nn=1:size(nc_files,1)

    nameFile=nc_files(nn).name;
    datevect{nn,1} = [nameFile(end-10:end-3),'000000'];   
    datevect_num(nn)=datenum([str2num(datevect{nn,1}(1:4)),str2num(datevect{nn,1}(5:6)),str2num(datevect{nn,1}(7:8)),0,0,0]);


    sea_ice_thickness = ncread(['./smos_sea_ice_thickness/' nameFile], 'sea_ice_thickness');
    vecice_frac(nn,1) = sea_ice_thickness(x,y);    

    uncertainty=ncread(['./smos_sea_ice_thickness/' nameFile], 'ice_thickness_uncertainty');
    ice_thick_uncertainty(nn)=uncertainty(x,y);


    Tsurf = ncread(['./smos_sea_ice_thickness/' nameFile], 'Tsurf');
    vecice_frac(nn,2) = Tsurf(x,y);
    if vecice_frac(nn,2) ==-999
        vecice_frac(nn,2) = NaN;
    else
        vecice_frac(nn,2) = vecice_frac(nn,2) -273.15; %%%% from Kelvin to celsius
    end


%     imagesc(sea_ice_thickness,'AlphaData',imAlpha);
%     set(gca,'color',[1 1 1]);
%     caxis([0 2])
%     colorbar
%     xlim([1 size(land,2)])
%     ylim([1 size(land,1)])
%     hold on
%     plot(y,x,'xk','markersize',16,'linewidth',3)
%     title(datestr([str2num(datevect{nn,1}(1:4)),str2num(datevect{nn,1}(5:6)),str2num(datevect{nn,1}(7:8)),0,0,0]))
%     pause(0.1)

end

T=table(cell2mat(datevect),vecice_frac(:,1),vecice_frac(:,2),'VariableNames',{'timestamp','sea_ice_thickness','Tsurf'});

writetable(T,[cd '/variables_SMOS_' sitename '.csv']) 

disp('Done')

%%
% 
% sea_ice_thickness(sea_ice_thickness<0)=NaN;
% 
% clc
% figure
% 
% %%% set axes
% latlimit=[65 90];
% centralmeridian=-160;
% maph=axesm('MapProjection','stereo','MapLatLimit',latlimit,  'Origin',[90 centralmeridian 0]);
% 
% %%% add grid
% parallel=[70 75 80 85];
% setm(maph,'Grid','on','Glinewidth',1,'PLineLocation',parallel)
% 
% %%% add grid labeling
% MLabelLocation=[-120 -60 0 60 120 180];
% Mpos=latlimit(1)-5;
% PLabelMeridian=30;
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
% land = shaperead('landareas.shp', 'UseGeoCoords', true);
% geoshow(maph, land, 'FaceColor',[0.80 0.80 0.80],'EdgeColor',0.30*[1 1 1]);
% %%% add ice data
% surfm(latitude,longitude,sea_ice_thickness)
% 
% 
% 
% %%
% 
% 
% 
