fileID = fopen([path_auxData 'north_x_y_lat_lon.txt']);
coord = fscanf(fileID, '%d %d %f %f',[4, Inf]);
fclose(fileID);
%%% Arctic data are 361*361
N=361;
latitude=reshape(coord(3,:,:),[N,N]);
longitude=reshape(coord(4,:,:),[N,N]);
%%%% map subset that has been chosen for the correlation computation
x_small=50:250;
y_small=100:250;
latitude_small=latitude(x_small, y_small);
longitude_small=longitude(x_small, y_small);

%%% reshape data
corr_small=NaN*zeros(length(x_small),length(y_small));
for ii=1:length(vecR)
    if vecP(ii)<0.05
        corr_small(vecR(ii,2),vecR(ii,3))=abs(vecR(ii,1)); 
    else
        corr_small(vecR(ii,2),vecR(ii,3))=0; 
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
surfm(latitude_small, longitude_small, corr_small)

%%% add mooring
gps_site = [72.90687 , -159.01807];
hold on
plotm(gps_site(1),gps_site(2),'xk','markersize',16,'linewidth',3)

colorbar