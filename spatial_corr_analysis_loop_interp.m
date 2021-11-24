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

Nvar_ok=size(auxData_ok, 2);

%% Limit to a single frequency band

load ANL_far

ff=2;
SPL_ANL_ok=SPL_ANL(:,ff);

Nosisaf=length(datenum_osisaf);
SPL_ANL_ave=NaN(size(datenum_osisaf));
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
        SPL_ANL_ave(tt)=mean(SPL_ANL_ok(i_min:i_max));
        SPL_ANL_ave2(tt)=mean(SPL_ANL_ok(i_min2:i_max2));
        SPL_ANL_ave2_far(tt)=nanmean(far_hor(i_min2:i_max2));
    end
end

figure,
plot(t,SPL_ANL_ok)
hold on
plot(datenum_osisaf, SPL_ANL_ave)
hold on
plot(datenum_osisaf,SPL_ANL_ave2)
hold on
plot(datenum_osisaf,SPL_ANL_ave2_far+20)
datetick


%% Interpolate
auxData_t_psd_near=interp1(datenum_osisaf, auxData_ok, t, 'nearest');
auxData_t_psd_linear=interp1(datenum_osisaf, auxData_ok, t, 'linear');

       
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


corr_spa_near=NaN([119 177 Nloop]);
pval_spa_near=NaN([119 177 Nloop]);
corr_spa_linear=NaN([119 177 Nloop]);
pval_spa_linear=NaN([119 177 Nloop]);
corr_spa_ave=NaN([119 177 Nloop]);
pval_spa_ave=NaN([119 177 Nloop]);
corr_spa_ave2=NaN([119 177 Nloop]);
pval_spa_ave2=NaN([119 177 Nloop]);

tt=1;
for t_num_loop=t0_num:15:t1_num
    t_beg_num=t_num_loop;
    t_end_num=t_num_loop+60;

    disp(datestr(t_beg_num, 'yyyymmdd'))
    %% Interpolated auxData

    ind_t_ok=find(t > t_beg_num & t < t_end_num);

%     t_ok=t(ind_t_ok);
%     Nt=length(t);
    auxData_t_psd_tt_near=auxData_t_psd_near(ind_t_ok,:);
    auxData_t_psd_tt_linear=auxData_t_psd_linear(ind_t_ok,:);
        
    vecSPL_tt=SPL_ANL_ok(ind_t_ok);

    %%%% Compute correlation

    for ii=1:Nvar_ok
        
        
        ind_no_nan=~isnan(auxData_t_psd_tt_near(:,ii));
        %%%% remove corr data which are nearly constant and/or nearly all
        %%%% NaN
        if length(unique(auxData_t_psd_tt_near(ind_no_nan,ii)))>15
            X = [vecSPL_tt(ind_no_nan), auxData_t_psd_tt_near(ind_no_nan,ii)];
            X_norma=zscore(X,[],1);
            [rho,pval] = corr(X_norma(:,1),X_norma(:,2),'type','Pearson','rows','all','tail','both');     
            
            ind_ini=ind_ice_ok(ii);
            [iii, jjj]=ind2sub(size(latitude), ind_ini);
            if pval < 0.05
                corr_spa_near(iii, jjj, tt)=rho;
            end  
            pval_spa_near(iii,jjj,tt)=pval;
        end
        J=squeeze(corr_spa_near(:,:,tt));
        [cmax_near(tt), indc_near(tt)]=max(J(:));
        
        ind_no_nan=~isnan(auxData_t_psd_tt_linear(:,ii));
        %%%% remove corr data which are nearly constant and/or nearly all
        %%%% NaN
        if length(unique(auxData_t_psd_tt_linear(ind_no_nan,ii)))>15
            X = [vecSPL_tt(ind_no_nan), auxData_t_psd_tt_linear(ind_no_nan,ii)];
            X_norma=zscore(X,[],1);
            [rho,pval] = corr(X_norma(:,1),X_norma(:,2),'type','Pearson','rows','all','tail','both');     
            
            ind_ini=ind_ice_ok(ii);
            [iii, jjj]=ind2sub(size(latitude), ind_ini);
            if pval < 0.05
                corr_spa_linear(iii, jjj, tt)=rho;
            end  
            pval_spa_linear(iii,jjj,tt)=pval;
        end
        J=squeeze(corr_spa_linear(:,:,tt));
        [cmax_linear(tt), indc_linear(tt)]=max(J(:));        
              
                    
    end

    
    
    
    %% Averaged SPL

    ind_t_ok=find(datenum_osisaf > t_beg_num & datenum_osisaf < t_end_num);

    auxData_tt=auxData_ok(ind_t_ok,:);      
    vecSPL_tt=SPL_ANL_ave(ind_t_ok);
    vecSPL_tt2=SPL_ANL_ave2(ind_t_ok);
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
            
            X2 = [vecSPL_tt2(ind_no_nan), auxData_tt(ind_no_nan,ii)];
            X_norma2=zscore(X2,[],1);
            [rho2,pval2] = corr(X_norma2(:,1),X_norma2(:,2),'type','Pearson','rows','all','tail','both');  
            if pval2 < 0.05
                corr_spa_ave2(iii, jjj, tt)=rho2;
            end  
            pval_spa_ave2(iii,jjj,tt)=pval2;           


        end
        J=squeeze(corr_spa_ave(:,:,tt));
        [cmax_ave(tt), indc_ave(tt)]=max(J(:));
              
        J2=squeeze(corr_spa_ave2(:,:,tt));
        [cmax_ave2(tt), indc_ave2(tt)]=max(J2(:));                                 
    end    
    

    %% Plot results
    figure('visible','off');
    h = get(0,'children');
    scrsz = get(0,'ScreenSize');
    set(h,'Position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)])  
    
%     subplot(221)   
%     %%% set axes
%     maph=axesm('MapProjection','lambertstd','MapLatLimit',latlimit,'MapLonLimit',lonlimit);
% 
%     %%% add grid
%     setm(maph,'Grid','on','Glinewidth',1,'PLineLocation',parallel, 'MLineLocation',MLabelLocation);
% 
%     %%% add grid labeling
%     setm(maph,'Fontangle','normal',...
%       'FontSize',12,'fontweight','b',...
%       'MeridianLabel','on',...
%       'MLabelLocation',MLabelLocation,...
%       'MLabelParallel',Mpos,...
%       'ParallelLabel','on',...
%       'PLabelLocation',parallel,...
%       'PLabelMeridian',PLabelMeridian);
% 
%     %%% add land
%     geoshow(maph, land, 'FaceColor',[0.80 0.80 0.80],'EdgeColor',0.30*[1 1 1]);
% 
%     %%% plot data
%     surfm(double(latitude), double(longitude), squeeze(corr_spa(:,:,tt)))
% 
%     %%% add mooring
%     plotm(gps_site(1),gps_site(2),'xk','markersize',16,'linewidth',3)
%     [toto, tata]=ind2sub(size(latitude), indc(tt));
%     plotm(latitude(toto, tata),longitude(toto,tata),'xr','markersize',16,'linewidth',3)
%     colorbar
%     caxis([0 0.5])
%     title(['Drift magnitude - Max corr = ' num2str(cmax(tt))])
    
    
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
    surfm(double(latitude), double(longitude), squeeze(corr_spa_near(:,:,tt)))

    %%% add mooring
    plotm(gps_site(1),gps_site(2),'xk','markersize',16,'linewidth',3)
    [toto, tata]=ind2sub(size(latitude), indc_near(tt));
    plotm(latitude(toto, tata),longitude(toto,tata),'xr','markersize',16,'linewidth',3)
    colorbar
    caxis([0 0.5])
    title(['Nearest interp - Max corr = ' num2str(cmax_near(tt))])
    
   
    
    
    
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
    surfm(double(latitude), double(longitude), squeeze(corr_spa_linear(:,:,tt)))

    %%% add mooring
    plotm(gps_site(1),gps_site(2),'xk','markersize',16,'linewidth',3)
    [toto, tata]=ind2sub(size(latitude), indc_linear(tt));
    plotm(latitude(toto, tata),longitude(toto,tata),'xr','markersize',16,'linewidth',3)
    colorbar
    caxis([0 0.5])
    title(['Linear interp - Max corr = ' num2str(cmax_linear(tt))])
    
    
    
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
    surfm(double(latitude), double(longitude), squeeze(corr_spa_ave(:,:,tt)))

    %%% add mooring
    plotm(gps_site(1),gps_site(2),'xk','markersize',16,'linewidth',3)
    [toto, tata]=ind2sub(size(latitude), indc_ave(tt));
    plotm(latitude(toto, tata),longitude(toto,tata),'xr','markersize',16,'linewidth',3)
    colorbar
    caxis([0 0.5])
    title(['1-day averaged SPL - Max corr = ' num2str(cmax_ave(tt))])   

    
    

    
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
    surfm(double(latitude), double(longitude), squeeze(corr_spa_ave2(:,:,tt)))

    %%% add mooring
    plotm(gps_site(1),gps_site(2),'xk','markersize',16,'linewidth',3)
    [toto, tata]=ind2sub(size(latitude), indc_ave(tt));
    plotm(latitude(toto, tata),longitude(toto,tata),'xr','markersize',16,'linewidth',3)
    colorbar
    caxis([0 0.5])
    title(['2-day averaged SPL - Max corr = ' num2str(cmax_ave2(tt))])    
    

    print(gcf,['./spatial_corr_result/' sitename '_interp/spatial_corr_' ...
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

save spatial_cor_results_interp.mat corr_spa_ave2 corr_spa_linear corr_spa_near corr_spa_ave date_loop f1 f2 latitude longitude SPL_ANL_ave2 SPL_ANL_ave SPL_ANL_ave2_far

disp('Done')