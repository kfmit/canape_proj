clear all
close all
clc

load ANL_SHRU5.mat
t=timestamp_num_spectro;
SPL_ANL_nobeamform=SPL_ANL;

load ANL_SHRU5_beamform.mat
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

%% Limit to a single frequency band 250-350 Hz
ff=5;
SPL_ANL_ok=SPL_ANL_nobeamform(:,ff);
ff=2;
SPL_ANL_up_ok=squeeze(SPL_ANL(:,ff,1));
SPL_ANL_hor_ok=squeeze(SPL_ANL(:,ff,3));

Nosisaf=length(datenum_osisaf);
SPL_ANL_ave=NaN(size(datenum_osisaf)); 
SPL_ANL_ave_up=NaN(size(datenum_osisaf)); 
SPL_ANL_ave_hor=NaN(size(datenum_osisaf)); 

for tt=1:Nosisaf-1
    [~, b1]=min(abs(t-datenum_osisaf(tt)));
    [~, b2]=min(abs(t-datenum_osisaf(tt+1)));
    db=b2-b1;
    if db>10
        i_min=b1-round(db);
        if i_min<1
            i_min=1;
        end
        i_max=min([length(t),b1+round(db)]);
        SPL_ANL_ave(tt)=mean(SPL_ANL_ok(i_min:i_max));
        SPL_ANL_ave_up(tt)=mean(SPL_ANL_up_ok(i_min:i_max));
        SPL_ANL_ave_hor(tt)=mean(SPL_ANL_hor_ok(i_min:i_max));
    end
end
% 
figure,
% plot(t,SPL_ANL_ok)
% hold on
plot(datenum_osisaf, SPL_ANL_ave)
hold on
plot(datenum_osisaf, SPL_ANL_ave_up)
hold on
plot(datenum_osisaf, SPL_ANL_ave_hor)
datetick
legend('Raw data', 'Beamforming up', 'Beamforming horizontally')
grid on

       
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


corr_spa_ave=NaN([119 177 Nloop]);
pval_spa_ave=NaN([119 177 Nloop]);
corr_spa_ave_up=NaN([119 177 Nloop]);
pval_spa_ave_up=NaN([119 177 Nloop]);
corr_spa_ave_hor=NaN([119 177 Nloop]);
pval_spa_ave_hor=NaN([119 177 Nloop]);

tt=1;
for t_num_loop=t0_num:15:t1_num
    t_beg_num=t_num_loop;
    t_end_num=t_num_loop+60;

    disp(datestr(t_beg_num, 'yyyymmdd'))



    ind_t_ok=find(datenum_osisaf > t_beg_num & datenum_osisaf < t_end_num);

    auxData_tt=auxData_ok(ind_t_ok,:);      
    vecSPL_tt=SPL_ANL_ave(ind_t_ok);
    vecSPL_tt_up=SPL_ANL_ave_up(ind_t_ok);
    vecSPL_tt_hor=SPL_ANL_ave_hor(ind_t_ok);
    %%%% Compute correlation

    for ii=1:Nvar_ok
     
        ind_no_nan=~isnan(auxData_tt(:,ii));
        %%%% remove corr data which are nearly constant and/or nearly all
        %%%% NaN
        if length(unique(auxData_tt(ind_no_nan,ii)))>10
            
            ind_ini=ind_ice_ok(ii);
            [iii, jjj]=ind2sub(size(latitude), ind_ini);
                        
          
            X = [vecSPL_tt(ind_no_nan), auxData_tt(ind_no_nan,ii)];
            X_norma=zscore(X,[],1);
            [rho,pval] = corr(X_norma(:,1),X_norma(:,2),'type','Pearson','rows','all','tail','both');    
            if pval < 0.05
                corr_spa_ave(iii, jjj, tt)=rho;
            end  
            pval_spa_ave(iii,jjj,tt)=pval;
            
            
            X = [vecSPL_tt_up(ind_no_nan), auxData_tt(ind_no_nan,ii)];
            X_norma=zscore(X,[],1);
            [rho,pval] = corr(X_norma(:,1),X_norma(:,2),'type','Pearson','rows','all','tail','both');    
            if pval < 0.05
                corr_spa_ave_up(iii, jjj, tt)=rho;
            end  
            pval_spa_ave_up(iii,jjj,tt)=pval;
            
            
            X = [vecSPL_tt_hor(ind_no_nan), auxData_tt(ind_no_nan,ii)];
            X_norma=zscore(X,[],1);
            [rho,pval] = corr(X_norma(:,1),X_norma(:,2),'type','Pearson','rows','all','tail','both');    
            if pval < 0.05
                corr_spa_ave_hor(iii, jjj, tt)=rho;
            end  
            pval_spa_ave_hor(iii,jjj,tt)=pval;            
            

        end
        J=squeeze(corr_spa_ave(:,:,tt));
        [cmax_ave(tt), indc_ave(tt)]=max(J(:));
        
        J=squeeze(corr_spa_ave_up(:,:,tt));
        [cmax_ave_up(tt), indc_ave_up(tt)]=max(J(:));
        
        J=squeeze(corr_spa_ave_hor(:,:,tt));
        [cmax_ave_hor(tt), indc_ave_hor(tt)]=max(J(:));
                               
    end    
    

    %% Plot results
    figure('visible','off');
    h = get(0,'children');
    scrsz = get(0,'ScreenSize');
    set(h,'Position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)])  
    

    
  
    
    
    subplot(221)   
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
    surfm(double(latitude), double(longitude), squeeze(corr_spa_ave(:,:,tt)))
    %%% add mooring
    plotm(gps_site(1),gps_site(2),'xk','markersize',16,'linewidth',3)
    [toto, tata]=ind2sub(size(latitude), indc_ave(tt));
    plotm(latitude(toto, tata),longitude(toto,tata),'xr','markersize',16,'linewidth',3)
    colorbar
    caxis([0 0.5])
    title(['No beamforming - Max corr = ' num2str(cmax_ave(tt))])   

    
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
    surfm(double(latitude), double(longitude), squeeze(corr_spa_ave_up(:,:,tt)))
    %%% add mooring
    plotm(gps_site(1),gps_site(2),'xk','markersize',16,'linewidth',3)
    [toto, tata]=ind2sub(size(latitude), indc_ave(tt));
    plotm(latitude(toto, tata),longitude(toto,tata),'xr','markersize',16,'linewidth',3)
    colorbar
    caxis([0 0.5])
    title(['Vertical beamforming - Max corr = ' num2str(cmax_ave_up(tt))])
    
    
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
    surfm(double(latitude), double(longitude), squeeze(corr_spa_ave_hor(:,:,tt)))
    %%% add mooring
    plotm(gps_site(1),gps_site(2),'xk','markersize',16,'linewidth',3)
    [toto, tata]=ind2sub(size(latitude), indc_ave(tt));
    plotm(latitude(toto, tata),longitude(toto,tata),'xr','markersize',16,'linewidth',3)
    colorbar
    caxis([0 0.5])
    title(['Horizontal beamforming - Max corr = ' num2str(cmax_ave_hor(tt))])    
    
    

    print(gcf,['./spatial_corr_result/' sitename '_beamform/spatial_corr_' ...
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

save spatial_cor_results_beamform.mat corr_spa_ave corr_spa_ave_up corr_spa_ave_hor date_loop f1 f2 latitude longitude SPL_ANL_ave SPL_ANL_ave_up SPL_ANL_ave_hor

disp('Done')