load lat_lon_osisaf


%%% reshape data
corr_spa=NaN(size(latitude));
for ii=1:length(vecR)
    if vecP(ii)<0.05
        corr_spa(vecR(ii,2),vecR(ii,3))=abs(vecR(ii,1)); 
    else
        corr_spa(vecR(ii,2),vecR(ii,3))=NaN; 
    end
end

figure('visible','off');


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
surfm(double(latitude), double(longitude), corr_spa)

%%% add mooring

if str2num(namefile(12))==1
    sitename = 'SHRU1';
    gps_site = [72+54.4123/60 , -(159+1.0840/60)];
elseif str2num(namefile(12))==2
    sitename = 'SHRU2';
    gps_site = [72+45.2347/60 , -(158+16.3243/60)];    
elseif str2num(namefile(12))==3
    sitename = 'SHRU3';
    gps_site = [72+40.6924/60 , -(157+54.6493/60)];
elseif str2num(namefile(12))==4
    sitename = 'SHRU4';
    gps_site = [72+36.6582/60 , -(157+32.2475/60)];    
elseif  str2num(namefile(12))==5
    sitename = 'SHRU5';
    gps_site = [72+54.4580/60 , -(157+29.2442/60)];  
else
    gps_site = [NaN, NaN];
end

plotm(gps_site(1),gps_site(2),'xk','markersize',16,'linewidth',3)

colorbar