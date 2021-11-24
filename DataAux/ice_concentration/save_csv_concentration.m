close all
clear all
clc

plot_map=0;
write_csv_data=1;


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


latitude = ncread(['./data_nc/' nc_files(1).name], 'latitude');
%%%% Satellite map has longitude between 360 and 0
longitude = ncread(['./data_nc/' nc_files(1).name], 'longitude');
longitude(longitude>180)=longitude(longitude>180)-360;


ind_ok=[];
err=0.05;
err_step=0.025;
while (isempty(ind_ok) && err<1)
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
land=ncread(['./data_nc/' nameFile], 'land');
ind_land=find(land);


latitude(x,y)
longitude(x,y)
figure
imagesc(land)
hold on
plot(y,x,'xk','markersize',16,'linewidth',3)

%%


imAlpha=ones(size(land));
imAlpha(ind_land)=0;

datevect={};
vecice_frac=zeros(size(nc_files,1),1);
for nn=1:size(nc_files,1)

    nameFile=nc_files(nn).name;
    ice_frac = ncread(['./data_nc/' nameFile], 'sea_ice_area_fraction');

    datevect{nn,1} = [nameFile(1:8),'000000'];
    datevect_num(nn)=datenum([str2num(datevect{nn,1}(1:4)),str2num(datevect{nn,1}(5:6)),str2num(datevect{nn,1}(7:8)),0,0,0]);
    vecice_frac(nn,1) = ice_frac(x,y);    
    
    if plot_map
        figure(1), hold on,
        imagesc(ice_frac,'AlphaData',imAlpha);
        set(gca,'color',[1 1 1]);
        caxis([0 100])
        colorbar
        xlim([1 size(land,2)])
        ylim([1 size(land,1)])
        plot(y,x,'xk','markersize',16,'linewidth',3)
        title(datestr([str2num(datevect{nn,1}(1:4)),str2num(datevect{nn,1}(5:6)),str2num(datevect{nn,1}(7:8)),0,0,0]))
        pause(0.05)
    end
end


%%
if write_csv_data
    T=table(cell2mat(datevect),vecice_frac,'VariableNames',[{'timestamp','icefrac'}]);
    writetable(T,[cd '/variables_ASI-SSMI_' sitename '.csv'])
end

    
%%
ice_frac=vecice_frac;


figure
plot(datevect_num, ice_frac)
grid on
ylabel('Ice fraction (%)')
datetick('x')
    
