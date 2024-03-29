%%%% A code to let me fuck around with the spatial correlation anlysis
% kat_spatial_corr_time_series_freqs and kat_spatial_corr_time_series_freq
% are the same thing


close all
clear all
clc

addpath('/home/kfung/Downloads/CANAPE/mat_files/')
addpath('/home/kfung/Downloads/CANAPE/spatial_mat_files/')
addpath('/home/kfung/Downloads/CANAPE/new_figs/')

%%% change name as neeeded %%%%%%%%%%%%%
% comes from spatial_corr_analysis_loop_interp_new_SHRU1.m
% ORIGINAL
% load spatial_cor_results_interp_new_shru1.mat
freq_range1 = 475;
freq_range2 = 525;
filename1 = ['spatial_cor_results_interp_new_shru1_50Hz_' num2str(freq_range1) '_' num2str(freq_range2) '.mat'];
filename2 = ['spatial_cor_results_interp_new_SHRU5_50Hz_' num2str(freq_range1) '_' num2str(freq_range2) '.mat'];

freq_range3 = 975;
freq_range4 = 1025;
filename3 = ['spatial_cor_results_interp_new_shru1_50Hz_' num2str(freq_range3) '_' num2str(freq_range4) '.mat'];
filename4 = ['spatial_cor_results_interp_new_SHRU5_50Hz_' num2str(freq_range3) '_' num2str(freq_range4) '.mat'];

% Load SHRU1 files, SKIPPED
%spatial_cor_results_interp_new_shru1_1250_1750.mat

% load(filename1);
% % save the specific variables to new names!
% corr_spa_ave2_shru1 = corr_spa_ave2;
% SPL_ANL_ave2_SHRU1 = SPL_ANL_ave2;
% gps_site = [72+54.4123/60 , -(159+1.0840/60)];
% 
% % SHRU 1
% load(filename3);
% % save the specific variables to new names!
% corr_spa_ave2_shru1_3 = corr_spa_ave2;
% SPL_ANL_ave2_SHRU1_3 = SPL_ANL_ave2;

%%%% comes from spatial_corr_analysis_loop_interp_SHRU5
%% LOAD SHRU5 Data
% load spatial_cor_results_interp_new.mat and spatial_cor_results_interp_new_1250_1750.mat

% NEW SHRU5 name
load(filename2);
% save the specific variables to new names!
corr_spa_ave2_shru5_2 = corr_spa_ave2;
SPL_ANL_ave2_SHRU5_2 = SPL_ANL_ave2;

load(filename4);
% save the specific variables to new names!
corr_spa_ave2_shru5_4 = corr_spa_ave2;
SPL_ANL_ave2_SHRU5_4 = SPL_ANL_ave2;



% GPS data
gps_site = [72+54.4580/60 , -(157+29.2442/60)];
gps_site_shru1 = [72+54.4123/60 , -(159+1.0840/60)];
dlon=40;
lonlimit=[gps_site(2)-dlon gps_site(2)+dlon];
lonlimit_ok=[lonlimit(2) lonlimit(1)+360];
latlimit=[65 85];


% comes from spatial_corr_analysis_loop_interp_new_SHRU5.m
load ANL_SHRU5_newfreq.mat % should already have everything
% load ANL_SHRU5.mat
t=timestamp_num_spectro;


% load auxData_icedrift_on_psd_time.mat
%%% contains lat/lon/dates/drift, doesn't need to change
% created from the icedrift csv code where?
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

% %% load ice edge and type
%
% ice_edge_dir='/home/julien/Desktop/DataAux/ice_edge/data_nc/ice_edge_nh_polstere-100_multi_';
% ddd='20161101';
% file_edge=[ice_edge_dir ddd '1200.nc'];
%
% latitude_edge = ncread(file_edge, 'lat');
% longitude_edge = ncread(file_edge, 'lon');
%
%
% ice_type_dir='/home/julien/Desktop/DataAux/ice_type/data_nc/ice_type_nh_polstere-100_multi_';
% ddd='20161101';
% file_edge=[ice_edge_dir ddd '1200.nc'];
%
% latitude_type = ncread(file_edge, 'lat');
% longitude_type = ncread(file_edge, 'lon');


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Loop for generating figure %%%%%%%%%%%%%%%%%%%%%%%
%%
% for tt=3:11, actual loop length is 11
loop_end = 11;
for tt=3:loop_end  % loop 1 long
    %     for tt=1
    t_beg_num=date_loop(1,tt);
    t_end_num=date_loop(2,tt);
    ind_t_ok_osisaf=find(t_osisaf >= t_beg_num & t_osisaf <= t_end_num);

    %% Read ice edge and type

    t_mid=(t_beg_num+t_end_num)/2;
    ddd_mid=(datestr(t_mid, 'yyyymmdd'));

    %     file_edge=[ice_edge_dir ddd_mid '1200.nc'];
    %     ice_edge = ncread(file_edge, 'ice_edge');
    %     edges=imgradient(ice_edge);
    %     edges(edges>.5)=10;
    %     toto_edge=find(edges==10);
    %
    %     latitude_edge_ok=double(latitude_edge(toto_edge));
    %     longitude_edge_ok=double(longitude_edge(toto_edge));
    %
    %
    %     file_type=[ice_type_dir ddd_mid '1200.nc'];
    %     ice_type = ncread(file_type, 'ice_type');
    %     types=imgradient(ice_type);
    %     types(types>.5)=10;
    %     toto_type=find(types==10);
    %
    %     latitude_type_ok=double(latitude_type(toto_type));
    %     longitude_type_ok=double(longitude_type(toto_type));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURE VISIBILITY HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure('visible','off');
    %%%%%%%%%%%%%%%%%%%%%%%%%%5%% FIGURE VISIBILITY HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     figure

    h = get(0,'children');
    scrsz = get(0,'ScreenSize');
%     set(h,'Position',[scrsz(1) scrsz(2) scrsz(3)/2 floor(scrsz(4)*0.66)]);
    set(h,'Position',scrsz);
%     set(h,'Position',[scrsz(1) scrsz(2) 2567 1112]); % [x y width height]
%     on full screen

    J=squeeze(corr_spa_ave2_shru5_2(:,:,tt));
    J(tutu)=NaN;
    %    figure, imagesc(J)
    [cmax(tt), indc_ave2(tt)]=max(abs(J(:)));
    [ii jj]=ind2sub(size(latitude), indc_ave2(tt));


    J_shru1=squeeze(corr_spa_ave2_shru5_2(:,:,tt));
    J_shru1(tutu)=NaN;

    %    figure, imagesc(J)
    [cmax_shru1(tt), indc_ave2_shru1(tt)]=max(abs(J_shru1(:)));
    [ii_shru1 jj_shru1]=ind2sub(size(latitude), indc_ave2_shru1(tt));


    %%%%%%%%%%%%%%%%%% SHRU5: BELOW  %%%%%%%%%%%%%%%%%%%
    %% Correlation map SHRU5, FREQ 1-2
    subplot(221)
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
    surfm(double(latitude), double(longitude), squeeze(corr_spa_ave2_shru5_2(:,:,tt)))

    %%% add mooring
    plotm(gps_site(1),gps_site(2),'xk','markersize',16,'linewidth',3)
    [toto, tata]=ind2sub(size(latitude), indc_ave2(tt));
    plotm(double(latitude(toto, tata)),double(longitude(toto,tata)),'xr','markersize',16,'linewidth',3)

    %     %%% add edges
    %     plotm(latitude_edge_ok, longitude_edge_ok, '.k','markersize',3,'linewidth',1)

    c=colorbar;
    c.Label.String = 'Correlation coefficient';
    caxis(ccc)
    %     title(['2-day averaged SPL - Max corr = ' num2str(R(tt))])
    title([datestr(t_beg_num, 'dd mmm yyyy') ' to ' datestr(t_end_num, 'dd mmm yyyy')], 'fontsize',20,'fontweight', 'bold')

    ylabel([num2str(freq_range1) '-' num2str(freq_range2)], 'fontsize',30,'fontweight', 'bold')

    dist(tt)=distance(gps_site(1), gps_site(2), double(latitude(toto, tata)),double(longitude(toto, tata)),referenceSphere('Earth'));

    %% Time series SHRU5, FREQ 1-2

    SPL_ANL_ok=SPL_ANL_ave2_SHRU5_2(ind_t_ok_osisaf);
    d_ok=squeeze(d(ind_t_ok_osisaf,ii,jj));
    %    d_ok_interp=interp1(t_osisaf(ind_t_ok_osisaf), d_ok, t(ind_t_ok), 'nearest');
    t_ok=t_osisaf(ind_t_ok_osisaf); %%% common time axis


    ind_no_nan=~isnan(d_ok);


    [SPL_ANL_ok_norm, mu_spl, sigma_spl]=zscore(SPL_ANL_ok(ind_no_nan));
    [d_ok_norm, mu_d, sigma_d]=zscore(d_ok(ind_no_nan));
    d_ok_norm2=(d_ok-mu_d)/sigma_d;

    %%%%%%%%% FINDING CORRELATION VALUES SHRU5 %%%%%%%%%%%%%%%%%
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

    subplot(222)
    plot(t_ok-t_ok(1),(SPL_ANL_ok-mu_spl)/sigma_spl, 'Color', [0    0.4470    0.7410], 'linewidth',2)

    hold on
    for bb=1:nb
        plot(t_ok(ind_b{bb})-t_ok(1), block{bb}, '-', 'Color', [0.8500    0.3250    0.0980], 'linewidth',2)
        hold on
    end

    plot(t_ok(ind_no_nan)-t_ok(1), d_ok_norm, 'o', 'Color', [0.8500    0.3250    0.0980], 'linewidth',2)
    title(['R_{max}=' num2str(R(tt))])
    xlabel('Days')
    grid on
    ylim([-3 4])
    legend(['ANL ' num2str(freq_range1) '-' num2str(freq_range2)], 'Ice drift')

    %%
    %%%%%%%%% Correlation map SHRU5, FREQ 3-5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(223)
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
    surfm(double(latitude), double(longitude), squeeze(corr_spa_ave2_shru5_4(:,:,tt)))

    %%% add mooring
    plotm(gps_site(1),gps_site(2),'xk','markersize',16,'linewidth',3)
    [toto, tata]=ind2sub(size(latitude), indc_ave2_shru1(tt));
    plotm(double(latitude(toto, tata)),double(longitude(toto,tata)),'xr','markersize',16,'linewidth',3)

    %     %%% add edges
    %     plotm(latitude_edge_ok, longitude_edge_ok, '.k','markersize',3,'linewidth',1)

    c=colorbar;
    c.Label.String = 'Correlation coefficient';
    caxis(ccc)
    %     title(['2-day averaged SPL - Max corr = ' num2str(R(tt))])
    %     title([datestr(t_beg_num, 'dd mmm yyyy') ' to ' datestr(t_end_num, 'dd mmm yyyy')], 'fontsize',20,'fontweight', 'bold')

    ylabel([num2str(freq_range3) '-' num2str(freq_range4)], 'fontsize',30,'fontweight', 'bold')

    dist_shru1(tt)=distance(gps_site_shru1(1), gps_site_shru1(2), double(latitude(toto, tata)),double(longitude(toto, tata)),referenceSphere('Earth'));

    %% Time series SHRU5, freq 3-4

    SPL_ANL_ok=SPL_ANL_ave2_SHRU5_4(ind_t_ok_osisaf);
    d_ok=squeeze(d(ind_t_ok_osisaf,ii_shru1,jj_shru1));
    %    d_ok_interp=interp1(t_osisaf(ind_t_ok_osisaf), d_ok, t(ind_t_ok), 'nearest');
    t_ok=t_osisaf(ind_t_ok_osisaf); %%% common time axis


    ind_no_nan=~isnan(d_ok);


    [SPL_ANL_ok_norm, mu_spl, sigma_spl]=zscore(SPL_ANL_ok(ind_no_nan));
    [d_ok_norm, mu_d, sigma_d]=zscore(d_ok(ind_no_nan));
    d_ok_norm2=(d_ok-mu_d)/sigma_d;
    
    %%%%%%%%%%% Correlation Values! %%%%%%%%%%%%%%%%
    [R_shru1(tt), Pvalue]=corr(SPL_ANL_ok_norm,d_ok_norm,'type','Pearson','rows','all','tail','both');

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



    subplot(224)
    plot(t_ok-t_ok(1),(SPL_ANL_ok-mu_spl)/sigma_spl, 'Color', [0    0.4470    0.7410], 'linewidth',2)

    hold on
    for bb=1:nb
        plot(t_ok(ind_b{bb})-t_ok(1), block{bb}, '-', 'Color', [0.8500    0.3250    0.0980], 'linewidth',2)
        hold on
    end

    plot(t_ok(ind_no_nan)-t_ok(1), d_ok_norm, 'o', 'Color', [0.8500    0.3250    0.0980], 'linewidth',2)
    title(['R_{max}=' num2str(R_shru1(tt))])
    xlabel('Days')
    grid on
    ylim([-3 4])

    legend(['ANL ' num2str(freq_range3) '-' num2str(freq_range4)], 'Ice drift')

    % title of the whole figure
    sgtitle(['Ice Drift Correlation for ' num2str(freq_range1) '-' num2str(freq_range2) ...
        ' to ' num2str(freq_range3) '-' num2str(freq_range4) ' Hz'])

    % printer currently OFF
    %%%%%%%%%%%% TURN ON AND OFF PRINTING %%%%%%%%%%%%%%%%%%%%%%%%
    % og path: './new_figs/spatial_corr_result/jasa_plot/'
        print(gcf,['./new_figs/spatial_corr_result_50/' num2str(freq_range1) '_' ... 
            num2str(freq_range4) '/spatial_corr_' datestr(t_beg_num, 'yyyymmdd') '-' datestr(t_end_num, 'yyyymmdd')]  ...
            ,'-dpng')

    %      print(gcf,['./spatial_corr_result/jasa_plot/spatial_corr_' ...
    %         datestr(t_beg_num, 'yyyymmdd') '-' datestr(t_end_num, 'yyyymmdd')]  ...
    %         ,'-dpng')

end
%%%% end of the figure generator %%%%


%%%%% Functions! %%%%%%%
dist_shrus=distance(gps_site(1), gps_site(2),gps_site_shru1(1), gps_site_shru1(2),referenceSphere('Earth'))/1000

dist_corr_shru5=[min(dist(3:loop_end)) mean(dist(3:loop_end)) max(dist(3:loop_end))]/1000

dist_corr_shru1=[min(dist_shru1(3:loop_end)) mean(dist_shru1(3:loop_end)) max(dist_shru1(3:loop_end))]/1000

r_shru5=[min(R(3:loop_end)) mean(R(3:loop_end)) max(R(3:loop_end))]

r_shru1=[min(R_shru1(3:loop_end)) mean(R_shru1(3:loop_end)) max(R_shru1(3:loop_end))]
