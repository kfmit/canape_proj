clear all
close all
clc

load ANL_SHRU5_newfreq.mat
t=timestamp_num_spectro;

addpath('/home/kfung/Downloads/CANAPE/auxData_SHRU5/')
addpath('/home/kfung/Downloads/CANAPE/DataAux/')

%% Geography limit
sitename = 'SHRU5';
gps_site = [72+54.4580/60 , -(157+29.2442/60)];  


latlimit=[65 85];
dlon=40;
lonlimit=[gps_site(2)-dlon gps_site(2)+dlon];
lonlimit_ok=[lonlimit(2) lonlimit(1)+360];
%% load env_data

load var_osisaf_new.mat
latitude=double(latitude);
longitude=double(longitude);
Nvar=119*177;

%%% restrict to area of interest
toto=find(latitude>=latlimit(1)+1 & latitude<=latlimit(2));
tata=find(longitude <= lonlimit_ok(1) | longitude >=lonlimit_ok(2));
ind_ice_ok=intersect(toto, tata);


auxData=reshape(d, [365,Nvar]);
auxData_ok=auxData(:,ind_ice_ok);

Nvar_ok=size(auxData_ok, 2);

%% load ice edge and type

ice_edge_dir='/home/kfung/Downloads/CANAPE/DataAux/ice_edge/data_nc/ice_edge_nh_polstere-100_multi_';
ddd='20161101';
file_edge=[ice_edge_dir ddd '1200.nc'];

latitude_edge = ncread(file_edge, 'lat');
longitude_edge = ncread(file_edge, 'lon');


ice_type_dir='/home/kfung/Downloads/CANAPE/DataAux/ice_type/data_nc/ice_type_nh_polstere-100_multi_';
ddd='20161101';
file_edge=[ice_edge_dir ddd '1200.nc'];

latitude_type = ncread(file_edge, 'lat');
longitude_type = ncread(file_edge, 'lon');

%% Limit to a single frequency band

load ANL_far

ff=4;
SPL_ANL_ok=SPL_ANL(:,ff);

Nosisaf=length(datenum_osisaf);
SPL_ANL_ave2=NaN(size(datenum_osisaf));
SPL_ANL_ave2_far=NaN(size(datenum_osisaf));

for tt=1:Nosisaf-1
    [~, b1]=min(abs(t-datenum_osisaf(tt)));
    [~, b2]=min(abs(t-datenum_osisaf(tt+1)));
    db=b2-b1;
    if db>10
        i_min=b1-round(db/2);
        i_min2=b1-round(db);
        if i_min<1
            i_min=1;
        end
        if i_min2<1
            i_min2=1;
        end
        i_max=min([length(t),b1+round(db/2)]);
        i_max2=min([length(t),b1+round(db)]);
        SPL_ANL_ave2(tt)=mean(SPL_ANL_ok(i_min2:i_max2));
        SPL_ANL_ave2_far(tt)=nanmean(far_hor(i_min2:i_max2));
    end
end

% figure,
% plot(t,SPL_ANL_ok)
% hold on
% plot(datenum_osisaf,SPL_ANL_ave2)
% hold on
% plot(datenum_osisaf,SPL_ANL_ave2_far)
% datetick

       
%% Prepare plots

%%%% prepare plot
latlimit=[65 85];
dlon=40;
lonlimit=[gps_site(2)-dlon gps_site(2)+dlon];
centralmeridian=-160;
parallel=[70 75 80];
MLabelLocation=[-180:20:180];
Mpos=latlimit(1)+2;
PLabelMeridian=30;
land = shaperead('landareas.shp', 'UseGeoCoords', true);
    
    
%% Must loop from Nov 1st to May 31

% t_beg=[2016,12,01,00,00,00];  %%% date vector [YY, MM, DD, h, min, sec]
% t_end=[2017,02,01,00,00,00];   %%% date vector [YY, MM, DD, h, min, sec]

t0=[2016,11,01,00,00,00];
t0_num=datenum(t0);

t1=[2017,05,31,00,00,00];
t1_num=datenum(t1);

Nloop=length(t0_num:15:t1_num);
date_loop=zeros(2,Nloop);

delta_loop=60;

corr_spa_ave2=NaN([119 177 Nloop]);
pval_spa_ave2=NaN([119 177 Nloop]);
corr_spa_ave2_far=NaN([119 177 Nloop]);
pval_spa_ave2_far=NaN([119 177 Nloop]);
corr_vec=NaN(1,Nvar_ok);
corr_vec_far=NaN(1,Nvar_ok);

tt=1;
for t_num_loop=t0_num:15:t1_num
    t_beg_num=t_num_loop;
    t_end_num=t_num_loop+delta_loop;

    ddd=(datestr(t_beg_num, 'yyyymmdd'));   
    disp(ddd)
    
    t_mid=(t_beg_num+t_end_num)/2;
    ddd_mid=(datestr(t_mid, 'yyyymmdd'));
    
    %% Read ice edge and type
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
    
    %% Averaged SPL

    ind_t_ok=find(datenum_osisaf > t_beg_num & datenum_osisaf < t_end_num);

    auxData_tt=auxData_ok(ind_t_ok,:);      
    vecSPL_tt2=SPL_ANL_ave2(ind_t_ok);
    vecSPL_tt2_far=SPL_ANL_ave2_far(ind_t_ok);
    %%%% Compute correlation

    for ii=1:Nvar_ok
     
        ind_no_nan=~isnan(auxData_tt(:,ii));
        %%%% remove corr data which are nearly constant and/or nearly all
        %%%% NaN
        if length(unique(auxData_tt(ind_no_nan,ii)))>15
            
            ind_ini=ind_ice_ok(ii);
            [iii, jjj]=ind2sub(size(latitude), ind_ini);
                        
            
            X2 = [vecSPL_tt2(ind_no_nan), auxData_tt(ind_no_nan,ii)];
            X_norma2=zscore(X2,[],1);
            [rho2,pval2] = corr(X_norma2(:,1),X_norma2(:,2),'type','Pearson','rows','all','tail','both');  
            if pval2 < 0.05
                corr_spa_ave2(iii, jjj, tt)=rho2;
            end  
            pval_spa_ave2(iii,jjj,tt)=pval2;   
            corr_vec(ii)=rho2;
            
            
            %%% there are a few Nan in vecSPL_tt2_far
            ind_no_nan_far=~isnan(auxData_tt(:,ii)) & ~isnan(vecSPL_tt2_far) ;
            X2_far = [vecSPL_tt2_far(ind_no_nan_far), auxData_tt(ind_no_nan_far,ii)];
            X_norma2_far=zscore(X2_far,[],1);
            [rho2_far,pval2_far] = corr(X_norma2_far(:,1),X_norma2_far(:,2),'type','Pearson','rows','all','tail','both');  
            if pval2_far < 0.05
                corr_spa_ave2_far(iii, jjj, tt)=rho2_far;
            end  
            pval_spa_ave2_far(iii,jjj,tt)=pval2_far; 
            corr_vec_far(ii)=rho2_far;

        end
              
        J2=squeeze(corr_spa_ave2(:,:,tt));
        [cmax_ave2(tt), indc_ave2(tt)]=max(J2(:));   
        [~, ind2]=max(abs(corr_vec));
                      
        J2_far=squeeze(corr_spa_ave2_far(:,:,tt));
        [cmax_ave2_far(tt), indc_ave2_far(tt)]=max(J2_far(:));  
        [~, ind2_far]=max(abs(corr_vec_far));
    end    
    

    %% Plot results
    figure('visible','off');
    h = get(0,'children');
    scrsz = get(0,'ScreenSize');
    set(h,'Position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)])  
    
    ccc=[0 0.7];  %%% for caxis
    
    subplot(221)
    plot((vecSPL_tt2-nanmean(vecSPL_tt2))/nanstd(vecSPL_tt2))
    hold on
    plot((auxData_tt(:,ind2)-nanmean(auxData_tt(:,ind2)))/nanstd(auxData_tt(:,ind2)))  
    grid on
    
    subplot(222)
    plot((vecSPL_tt2_far-nanmean(vecSPL_tt2_far))/nanstd(vecSPL_tt2_far))
    hold on
    plot((auxData_tt(:,ind2)-nanmean(auxData_tt(:,ind2)))/nanstd(auxData_tt(:,ind2)))  
    grid on    
    
    
    
    subplot(223)   
    %%% set axes
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
    plotm(latitude(toto, tata),longitude(toto,tata),'xr','markersize',16,'linewidth',3)
    colorbar
    caxis(ccc)
    title(['2-day averaged SPL - Max corr = ' num2str(cmax_ave2(tt))])    
    %%% add edges
    plotm(latitude_edge_ok, longitude_edge_ok, '.k','markersize',5,'linewidth',1)
    %%% add type
%     plotm(latitude_type_ok, longitude_type_ok, '.r','markersize',5,'linewidth',1)    
    
    
    subplot(224)   
    %%% set axes
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
    surfm(double(latitude), double(longitude), squeeze(corr_spa_ave2_far(:,:,tt)))

    %%% add mooring
    plotm(gps_site(1),gps_site(2),'xk','markersize',16,'linewidth',3)
    [toto, tata]=ind2sub(size(latitude), indc_ave2_far(tt));
    plotm(latitude(toto, tata),longitude(toto,tata),'xr','markersize',16,'linewidth',3)
    colorbar
    caxis(ccc)
    title(['2-day averaged SPL FAR - Max corr = ' num2str(cmax_ave2_far(tt))])      
    %%% add edges
    plotm(latitude_edge_ok, longitude_edge_ok, '.k','markersize',5,'linewidth',1)
    %%% add type
%     plotm(latitude_type_ok, longitude_type_ok, '.r','markersize',5,'linewidth',1)    
    
    
    print(gcf,['./new_figs/spatial_corr_result/' sitename '_interp_new/spatial_corr_' ...
        '-freq_' num2str(f1(ff)) '-' num2str(f2(ff)) ...
        '_' datestr(t_beg_num, 'yyyymmdd') '-' datestr(t_end_num, 'yyyymmdd')]  ...
        ,'-dpng')

%         saveas(gcf,['./spatial_corr_result/' sitename '/spatial_corr_' ...
%             '-freq_' num2str(f1(ff)) '-' num2str(f2(ff)) ...
%             '_' datestr(t_beg_num, 'yyyymmdd') '-' datestr(t_end_num, 'yyyymmdd')],'fig')


    date_loop(1,tt)=t_beg_num;
    date_loop(2,tt)=t_end_num;
    tt=tt+1;  
    
    
end

save spatial_cor_results_interp_new_1250_1750.mat corr_spa_ave2 corr_spa_ave2_far date_loop f1 f2 latitude longitude SPL_ANL_ave2 SPL_ANL_ave2_far

disp('Done')