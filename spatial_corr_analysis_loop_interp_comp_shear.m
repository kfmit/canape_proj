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
    
  
    %% Averaged SPL

    ind_t_ok=find(datenum_osisaf > t_beg_num & datenum_osisaf < t_end_num);

    auxData_tt=auxData_ok(ind_t_ok,:);      
    vecSPL_tt2=SPL_ANL_ave2(ind_t_ok);
    vecSPL_tt2_far=SPL_ANL_ave2_far(ind_t_ok);
    
%     d_shear_tt=squeeze(nanmean(d_shear(ind_t_ok,:,:),1));
%     d_curl_tt=squeeze(nanmean(d_curl(ind_t_ok,:,:),1));
%     d_div_tt=squeeze(nanmean(d_div(ind_t_ok,:,:),1));

    d_shear_tt=squeeze(d_shear(ind_t_ok(end/2),:,:));
    d_curl_tt=squeeze(d_curl(ind_t_ok(end/2),:,:));
    d_div_tt=squeeze(d_div(ind_t_ok(end/2),:,:));
    
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
            
            

        end
              
        J2=squeeze(corr_spa_ave2(:,:,tt));
        [cmax_ave2(tt), indc_ave2(tt)]=max(J2(:));   
        [~, ind2]=max(abs(corr_vec));
                      
    end    
    

    %% Plot results
    figure('visible','off');
    h = get(0,'children');
    scrsz = get(0,'ScreenSize');
    set(h,'Position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)])  
    ccc=[0 0.7];  %%% for caxis for correlation
    ccc2=[0 3];  %%% for caxis div curl shear
   
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
    surfm(double(latitude), double(longitude), squeeze(corr_spa_ave2(:,:,tt)))

    %%% add mooring
    plotm(gps_site(1),gps_site(2),'xk','markersize',16,'linewidth',3)
    [toto, tata]=ind2sub(size(latitude), indc_ave2(tt));
    plotm(latitude(toto, tata),longitude(toto,tata),'xr','markersize',16,'linewidth',3)
    colorbar
    caxis(ccc)
    title(['2-day averaged SPL - Max corr = ' num2str(cmax_ave2(tt))])    
  
    
   
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
    surfm(double(latitude), double(longitude), abs(d_shear_tt))

    %%% add mooring
    plotm(gps_site(1),gps_site(2),'xk','markersize',16,'linewidth',3)
    [toto, tata]=ind2sub(size(latitude), indc_ave2(tt));
    plotm(latitude(toto, tata),longitude(toto,tata),'xr','markersize',16,'linewidth',3)
    colorbar
    caxis(ccc2)
    title('Mean shear of the drift')   

    
    
    
    
    
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
    surfm(double(latitude), double(longitude), abs(d_div_tt))

    %%% add mooring
    plotm(gps_site(1),gps_site(2),'xk','markersize',16,'linewidth',3)
    [toto, tata]=ind2sub(size(latitude), indc_ave2(tt));
    plotm(latitude(toto, tata),longitude(toto,tata),'xr','markersize',16,'linewidth',3)
    colorbar
    caxis(ccc2)
    title('Mean div of the drift') 
    
    
    
   
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
    surfm(double(latitude), double(longitude), abs(d_curl_tt))

    %%% add mooring
    plotm(gps_site(1),gps_site(2),'xk','markersize',16,'linewidth',3)
    [toto, tata]=ind2sub(size(latitude), indc_ave2(tt));
    plotm(latitude(toto, tata),longitude(toto,tata),'xr','markersize',16,'linewidth',3)
    colorbar
    caxis(ccc2)
    title('Mean curl of the drift')   
    
     
    
    print(gcf,['./spatial_corr_result/' sitename '_interp_comp_shear/spatial_corr_' ...
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


disp('Done')