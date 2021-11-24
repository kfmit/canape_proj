close all
clear all
clc


load spatial_cor_results_interp_new.mat
gps_site = [72+54.4580/60 , -(157+29.2442/60)];  
dlon=40;
lonlimit=[gps_site(2)-dlon gps_site(2)+dlon];
lonlimit_ok=[lonlimit(2) lonlimit(1)+360];
latlimit=[65 85];

load ANL_SHRU5.mat
t=timestamp_num_spectro;

% load auxData_icedrift_on_psd_time.mat

load var_osisaf.mat
t_osisaf=datenum_osisaf;

% %% Reject area of no-interest
toto=find(latitude<=latlimit(1) | latitude>=latlimit(2));
tata=find(longitude >= lonlimit_ok(1) & longitude <=lonlimit_ok(2));
tutu=union(toto, tata);
% %%%%% verification
% % lat_ok=latitude;
% % lon_ok=longitude;
% % lat_ok(tutu)=NaN;
% % lon_ok(tutu)=NaN;
% % figure, imagesc(lat_ok)
% % figure, imagesc(lon_ok)


Nloop=size(date_loop,2);

%% load ice edge and type

ice_edge_dir='/home/julien/Desktop/DataAux/ice_edge/data_nc/ice_edge_nh_polstere-100_multi_';
ddd='20161101';
file_edge=[ice_edge_dir ddd '1200.nc'];

latitude_edge = ncread(file_edge, 'lat');
longitude_edge = ncread(file_edge, 'lon');


ice_type_dir='/home/julien/Desktop/DataAux/ice_type/data_nc/ice_type_nh_polstere-100_multi_';
ddd='20161101';
file_edge=[ice_edge_dir ddd '1200.nc'];

latitude_type = ncread(file_edge, 'lat');
longitude_type = ncread(file_edge, 'lon');


%% prepare plot
latlimit=[65 85];
dlon=40;
lonlimit=[gps_site(2)-dlon gps_site(2)+dlon];
centralmeridian=-160;
parallel=[70 75 80];
MLabelLocation=[-180:20:180];
Mpos=latlimit(1)+2;
PLabelMeridian=30;
land = shaperead('landareas.shp', 'UseGeoCoords', true);

ccc=[0 0.7];  %%% for caxis

%% Loop

for tt=1:Nloop
   
    
   
   
   t_beg_num=date_loop(1,tt);
   t_end_num=date_loop(2,tt);      
   ind_t_ok_osisaf=find(t_osisaf >= t_beg_num & t_osisaf <= t_end_num);
    
   
   %% Read ice edge and type
    
    t_mid=(t_beg_num+t_end_num)/2;
    ddd_mid=(datestr(t_mid, 'yyyymmdd'));
    
    file_edge=[ice_edge_dir ddd_mid '1200.nc'];
    ice_edge = ncread(file_edge, 'ice_edge');
    edges=imgradient(ice_edge);
    edges(edges>.5)=10;
    toto_edge=find(edges==10);
    
    latitude_edge_ok=double(latitude_edge(toto_edge));
    longitude_edge_ok=double(longitude_edge(toto_edge));
    
    
    file_type=[ice_type_dir ddd_mid '1200.nc'];
    ice_type = ncread(file_type, 'ice_type');
    types=imgradient(ice_type);
    types(types>.5)=10;
    toto_type=find(types==10);
    
    latitude_type_ok=double(latitude_type(toto_type));
    longitude_type_ok=double(longitude_type(toto_type));  
    
    
    
    
    
    figure('visible','off');
    h = get(0,'children');
    scrsz = get(0,'ScreenSize');
    set(h,'Position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)/2]) 
    
   J=squeeze(corr_spa_ave2(:,:,tt));
   J(tutu)=NaN;
%    figure, imagesc(J)  
   [cmax(tt), indc_ave2(tt)]=max(abs(J(:)));
   [ii jj]=ind2sub(size(latitude), indc_ave2(tt));
   
   
   
     %% Correlation map
    subplot(121)   
    maph=axesm('MapProjection','lambertstd','MapLatLimit',latlimit,'MapLonLimit',lonlimit);

    %%% add grid
    setm(maph,'Grid','on','Glinewidth',1,'PLineLocation',parallel, 'MLineLocation',MLabelLocation);

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

    %%% plot data
    surfm(double(latitude), double(longitude), squeeze(corr_spa_ave2(:,:,tt)))

    %%% add mooring
    plotm(gps_site(1),gps_site(2),'xk','markersize',16,'linewidth',3)
    [toto, tata]=ind2sub(size(latitude), indc_ave2(tt));
    plotm(double(latitude(toto, tata)),double(longitude(toto,tata)),'xr','markersize',16,'linewidth',3)
    
    %%% add edges
    plotm(latitude_edge_ok, longitude_edge_ok, '.k','markersize',5,'linewidth',1)
    
    colorbar
    caxis(ccc)
%     title(['2-day averaged SPL - Max corr = ' num2str(R(tt))])  
    title([datestr(t_beg_num, 'dd mmm yyyy') ' to ' datestr(t_end_num, 'dd mmm yyyy')])
    

    dist(tt)=distance(gps_site(1), gps_site(2), double(latitude(toto, tata)),double(longitude(toto, tata)),referenceSphere('Earth'));
   
   %% Time series   
   
     
  
   SPL_ANL_ok=SPL_ANL_ave2(ind_t_ok_osisaf);
   d_ok=squeeze(d(ind_t_ok_osisaf,ii,jj));   
%    d_ok_interp=interp1(t_osisaf(ind_t_ok_osisaf), d_ok, t(ind_t_ok), 'nearest');
   t_ok=t_osisaf(ind_t_ok_osisaf); %%% common time axis
   
      
    ind_no_nan=~isnan(d_ok);
   

    [SPL_ANL_ok_norm, mu_spl, sigma_spl]=zscore(SPL_ANL_ok(ind_no_nan));
    [d_ok_norm, mu_d, sigma_d]=zscore(d_ok(ind_no_nan));
    d_ok_norm2=(d_ok-mu_d)/sigma_d;
    
    [R(tt), Pvalue]=corr(SPL_ANL_ok_norm,d_ok_norm,'type','Pearson','rows','all','tail','both');
    
    toto=find(~isnan(d_ok));
    
    Ntoto=length(toto);
    nb=1;
    block{nb}=d_ok_norm2(toto(1));
    ind_b{nb}=toto(1);
    for iii=2:Ntoto
        if toto(iii)-toto(iii-1) == 1
            block{nb}=[block{nb} d_ok_norm2(toto(iii))];
            ind_b{nb}=[ind_b{nb} toto(iii)];
        else %%% new block
            nb=nb+1;
            block{nb}= d_ok_norm2(toto(iii));     
            ind_b{nb}= toto(iii);
        end
    end
                
    

    subplot(122)
    plot(t_ok-t_ok(1),(SPL_ANL_ok-mu_spl)/sigma_spl, 'Color', [0    0.4470    0.7410])
       
    hold on
    for bb=1:nb
        plot(t_ok(ind_b{bb})-t_ok(1), block{bb}, '-', 'Color', [0.8500    0.3250    0.0980])
        hold on
    end
    
    plot(t_ok(ind_no_nan)-t_ok(1), d_ok_norm, 'o', 'Color', [0.8500    0.3250    0.0980])
    title(['R_{max}=' num2str(R(tt))])
    xlabel('Days')
    grid on
    
        print(gcf,['./spatial_corr_result/plot_serenade/spatial_corr_' ...
        '_' datestr(t_beg_num, 'yyyymmdd') '-' datestr(t_end_num, 'yyyymmdd')]  ...
        ,'-dpng')
end

  
