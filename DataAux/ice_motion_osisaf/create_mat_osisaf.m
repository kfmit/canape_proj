close all
clear all
clc

plot_map=0;
write_csv_data=1;



nc_files = dir('./data/**/*.nc');


nn=105;
nameFile=nc_files(nn).name;
pathFile=nc_files(nn).folder;

latitude=ncread([pathFile '/' nameFile],'lat');
longitude=ncread([pathFile '/' nameFile],'lon');
%%%% Ice data has longitude between -180 and 180


%% Loop
Ndate=size(nc_files,1);
d=NaN(Ndate,119,177);
d_div=NaN(Ndate,119,177);
d_curl=NaN(Ndate,119,177);
d_shear=NaN(Ndate,119,177);

date_num_ok=zeros(Ndate,1);

datevect={};
for nn=1:Ndate
    nameFile=nc_files(nn).name;
    pathFile=nc_files(nn).folder;

    date_vec=nameFile(end-14:end-3);
    date_num=datenum(date_vec,'yyyymmddHHMM');
    date_num_ok(nn)=date_num-1;
    date_vec_ok=datestr(date_num_ok(nn),'yyyymmddHHMMSS');
    datevect{nn,1} = date_vec_ok;
    
    dx=ncread([pathFile '/' nameFile],'dX');
    dy=ncread([pathFile '/' nameFile],'dY');
    d(nn,:,:)=sqrt(dx.^2+dy.^2);
    d_div(nn,:,:)=divergence(dx,dy);
    [d_curl(nn,:,:) d_shear(nn,:,:)]=shear_2d(dx,dy);
    
    
    if plot_map
        figure(1)
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
        surfm(double(latitude), double(longitude), squeeze(d(nn,:,:)))
        title(date_vec_ok)
        pause(0.01)
    end
            
    
end

datenum_osisaf=date_num_ok;

save var_osisaf_new d d_div d_curl d_shear latitude longitude datenum_osisaf
