clear all
close all
clc

load ANL_SHRU5.mat
t=timestamp_num_spectro;

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

auxData_curl=reshape(d_curl, [365,Nvar]);
auxData_ok_curl=auxData_curl(:,ind_ice_ok);

auxData_div=reshape(d_div, [365,Nvar]);
auxData_ok_div=auxData_div(:,ind_ice_ok);

auxData_shear=reshape(d_shear, [365,Nvar]);
auxData_ok_shear=auxData_shear(:,ind_ice_ok);

Nvar_ok=size(auxData_ok, 2);

%% Limit to a single frequency band
ff=2;
SPL_ANL_ok=SPL_ANL(:,ff);

%% Interpolate
auxData_t_psd=interp1(datenum_osisaf, auxData_ok, t);
auxData_t_psd_curl=interp1(datenum_osisaf, auxData_ok_curl, t);
auxData_t_psd_div=interp1(datenum_osisaf, auxData_ok_div, t);
auxData_t_psd_shear=interp1(datenum_osisaf, auxData_ok_shear, t);
       
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

corr_spa=NaN([119 177 Nloop]);
pval_spa=NaN([119 177 Nloop]);
corr_spa_curl=NaN([119 177 Nloop]);
pval_spa_curl=NaN([119 177 Nloop]);
corr_spa_div=NaN([119 177 Nloop]);
pval_spa_div=NaN([119 177 Nloop]);
corr_spa_shear=NaN([119 177 Nloop]);
pval_spa_shear=NaN([119 177 Nloop]);

tt=1;
for t_num_loop=t0_num:15:t1_num
    t_beg_num=t_num_loop;
    t_end_num=t_num_loop+60;

    disp(datestr(t_beg_num, 'yyyymmdd'))
    %% Restric to time of interest


    ind_t_ok=find(t > t_beg_num & t < t_end_num);

    t_ok=t(ind_t_ok);
    Nt=length(t);
    auxData_t_psd_tt=auxData_t_psd(ind_t_ok,:);
    auxData_t_psd_tt_curl=auxData_t_psd_curl(ind_t_ok,:);
    auxData_t_psd_tt_div=auxData_t_psd_div(ind_t_ok,:);
    auxData_t_psd_tt_shear=auxData_t_psd_shear(ind_t_ok,:);
    
    
    vecSPL_tt=SPL_ANL_ok(ind_t_ok);

    %% Compute correlation

    for ii=1:Nvar_ok
        ind_no_nan=~isnan(auxData_t_psd_tt(:,ii));
        %%%% remove corr data which are nearly constant and/or nearly all
        %%%% NaN
        if length(unique(auxData_t_psd_tt(ind_no_nan,ii)))>15
            X = [vecSPL_tt(ind_no_nan), auxData_t_psd_tt(ind_no_nan,ii)];
            X_norma=zscore(X,[],1);
            [rho,pval] = corr(X_norma(:,1),X_norma(:,2),'type','Pearson','rows','all','tail','both');     
            
            ind_ini=ind_ice_ok(ii);
            [iii, jjj]=ind2sub(size(latitude), ind_ini);
            if pval < 0.05
                corr_spa(iii, jjj, tt)=rho;
            end  
            pval_spa(iii,jjj,tt)=pval;
        end
        J=squeeze(corr_spa(:,:,tt));
        [cmax(tt), indc(tt)]=max(J(:));
        
        
        ind_no_nan=~isnan(auxData_t_psd_tt_curl(:,ii));
        %%%% remove corr data which are nearly constant and/or nearly all
        %%%% NaN
        if length(unique(auxData_t_psd_tt_curl(ind_no_nan,ii)))>15
            X = [vecSPL_tt(ind_no_nan), auxData_t_psd_tt_curl(ind_no_nan,ii)];
            X_norma=zscore(X,[],1);
            [rho,pval] = corr(X_norma(:,1),X_norma(:,2),'type','Pearson','rows','all','tail','both');     
            
            ind_ini=ind_ice_ok(ii);
            [iii, jjj]=ind2sub(size(latitude), ind_ini);
            if pval < 0.05
                corr_spa_curl(iii, jjj, tt)=rho;
            end  
            pval_spa_curl(iii,jjj,tt)=pval;
        end
        J=squeeze(corr_spa_curl(:,:,tt));
        [cmax_curl(tt), indc_curl(tt)]=max(J(:));
        
        ind_no_nan=~isnan(auxData_t_psd_tt_div(:,ii));
        %%%% remove corr data which are nearly constant and/or nearly all
        %%%% NaN
        if length(unique(auxData_t_psd_tt_div(ind_no_nan,ii)))>15
            X = [vecSPL_tt(ind_no_nan), auxData_t_psd_tt_div(ind_no_nan,ii)];
            X_norma=zscore(X,[],1);
            [rho,pval] = corr(X_norma(:,1),X_norma(:,2),'type','Pearson','rows','all','tail','both');     
            
            ind_ini=ind_ice_ok(ii);
            [iii, jjj]=ind2sub(size(latitude), ind_ini);
            if pval < 0.05
                corr_spa_div(iii, jjj, tt)=rho;
            end  
            pval_spa_div(iii,jjj,tt)=pval;
        end
        J=squeeze(corr_spa_div(:,:,tt));
        [cmax_div(tt), indc_div(tt)]=max(J(:));        
        
        ind_no_nan=~isnan(auxData_t_psd_tt_shear(:,ii));
        %%%% remove corr data which are nearly constant and/or nearly all
        %%%% NaN
        if length(unique(auxData_t_psd_tt_shear(ind_no_nan,ii)))>15
            X = [vecSPL_tt(ind_no_nan), auxData_t_psd_tt_shear(ind_no_nan,ii)];
            X_norma=zscore(X,[],1);
            [rho,pval] = corr(X_norma(:,1),X_norma(:,2),'type','Pearson','rows','all','tail','both');     
            
            ind_ini=ind_ice_ok(ii);
            [iii, jjj]=ind2sub(size(latitude), ind_ini);
            if pval < 0.05
                corr_spa_shear(iii, jjj, tt)=rho;
            end  
            pval_spa_shear(iii,jjj,tt)=pval;
        end
        J=squeeze(corr_spa_div(:,:,tt));
        [cmax_shear(tt), indc_shear(tt)]=max(J(:));         
              
        
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
    surfm(double(latitude), double(longitude), squeeze(corr_spa(:,:,tt)))

    %%% add mooring
    plotm(gps_site(1),gps_site(2),'xk','markersize',16,'linewidth',3)
    [toto, tata]=ind2sub(size(latitude), indc(tt));
    plotm(latitude(toto, tata),longitude(toto,tata),'xr','markersize',16,'linewidth',3)
    colorbar
    caxis([0 0.5])
    title(['Drift magnitude - Max corr = ' num2str(cmax(tt))])
    
    
    subplot(222)   
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
    surfm(double(latitude), double(longitude), squeeze(corr_spa_curl(:,:,tt)))

    %%% add mooring
    plotm(gps_site(1),gps_site(2),'xk','markersize',16,'linewidth',3)
    [toto, tata]=ind2sub(size(latitude), indc_curl(tt));
    plotm(latitude(toto, tata),longitude(toto,tata),'xr','markersize',16,'linewidth',3)
    colorbar
    caxis([0 0.5])
    title(['Drift curl - Max corr = ' num2str(cmax_curl(tt))])
    
   
    
    
    
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
    surfm(double(latitude), double(longitude), squeeze(corr_spa_div(:,:,tt)))

    %%% add mooring
    plotm(gps_site(1),gps_site(2),'xk','markersize',16,'linewidth',3)
    [toto, tata]=ind2sub(size(latitude), indc_div(tt));
    plotm(latitude(toto, tata),longitude(toto,tata),'xr','markersize',16,'linewidth',3)
    colorbar
    caxis([0 0.5])
    title(['Drift divergence - Max corr = ' num2str(cmax_div(tt))])
    
    
    
    
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
    surfm(double(latitude), double(longitude), squeeze(corr_spa_shear(:,:,tt)))

    %%% add mooring
    plotm(gps_site(1),gps_site(2),'xk','markersize',16,'linewidth',3)
    [toto, tata]=ind2sub(size(latitude), indc_shear(tt));
    plotm(latitude(toto, tata),longitude(toto,tata),'xr','markersize',16,'linewidth',3)
    colorbar
    caxis([0 0.5])
    title(['Drift divergence - Max corr = ' num2str(cmax_shear(tt))])
    

    print(gcf,['./spatial_corr_result/' sitename '_shear/spatial_corr_' ...
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

save spatial_cor_results_shear.mat corr_spa corr_spa_shear corr_spa_div corr_spa_curl date_loop f1 f2 latitude longitude 

disp('Done')