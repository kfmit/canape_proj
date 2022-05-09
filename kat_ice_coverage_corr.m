%%%% A code to let me fuck around with the spatial correlation anlysis
% For correlation maps
% Make a threshold and look at all the pixels with R>5
% Report number of pixels (more or less the size of the area of interest)
% Report min, median and average distance between the pixels and the hydrophone 



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

freq_array1=[40 250 450 900 1250];
freq_array2=[60 350 550 1100 1750];

for i = 1:5 % length(freq_array2)

freq_range1 = freq_array1(i);
freq_range2 = freq_array2(i);
filename1 = ['spatial_cor_results_interp_new_shru1_' num2str(freq_range1) '_' num2str(freq_range2) '.mat'];
filename2 = ['spatial_cor_results_interp_new_' num2str(freq_range1) '_' num2str(freq_range2) '.mat'];
% 
freq_range3 = 1250;
freq_range4 = 1750;
filename3 = ['spatial_cor_results_interp_new_shru1_' num2str(freq_range3) '_' num2str(freq_range4) '.mat'];
filename4 = ['spatial_cor_results_interp_new_' num2str(freq_range3) '_' num2str(freq_range4) '.mat'];
% NEW
load(filename1);
%spatial_cor_results_interp_new_shru1_1250_1750.mat

corr_spa_ave2_shru1 = corr_spa_ave2;
SPL_ANL_ave2_SHRU1 = SPL_ANL_ave2;
gps_site_shru1 = [72+54.4123/60 , -(159+1.0840/60)];

%%%% comes from spatial_corr_analysis_loop_interp_SHRU5
% ORIGINAL
% load spatial_cor_results_interp_new.mat

% NEW
load(filename2);
%spatial_cor_results_interp_new_1250_1750.mat

gps_site = [72+54.4580/60 , -(157+29.2442/60)];
gps_site_shru1 = [72+54.4123/60 , -(159+1.0840/60)];
dlon=40;
lonlimit=[gps_site(2)-dlon gps_site(2)+dlon];
lonlimit_ok=[lonlimit(2) lonlimit(1)+360];
latlimit=[65 85];


% comes from spatial_corr_analysis_loop_interp_new_SHRU5.m
load ANL_SHRU5_newfreq.mat
% load ANL_SHRU5.mat
t=timestamp_num_spectro;
%%% WHERE DO I CONTROL FREQ SPECTRUM

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

% change labelling of map here
latlimit=[65 85];
dlon=40;
lonlimit=[gps_site(2)-dlon gps_site(2)+dlon];
centralmeridian=-160;
parallel=[75 80];
MLabelLocation=[-180:30:180];
Mpos=latlimit(1)+2;
PLabelMeridian=30;
land = shaperead('landareas.shp', 'UseGeoCoords', true);

ccc=[0 0.7];  %%% for caxis, will have to modiy to match the 0.5 cutoff AFTER

%% Loop

loop_end = 11;
for tt=3:loop_end
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
%     tiledlayout(1,2)
    %%%%%%%%%%%%%%%%%%%%%%%%%%5%% FIGURE VISIBILITY HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     figure
%     h = get(0,'children');
%     scrsz = get(0,'ScreenSize');
%     set(h,'Position',[scrsz(1) scrsz(2) scrsz(3)/2 floor(scrsz(4)*0.66)])
    
    % squeeze removes the dimension length 1, J become lat & lon for time
    % tt
    J=squeeze(corr_spa_ave2(:,:,tt));
    J(tutu)=NaN;

   

    %    figure, imagesc(J)
    [cmax(tt), indc_ave2(tt)]=max(abs(J(:)));
    [ii jj]=ind2sub(size(latitude), indc_ave2(tt));


    J_shru1=squeeze(corr_spa_ave2_shru1(:,:,tt));
    J_shru1(tutu)=NaN;

    %    figure, imagesc(J)
    [cmax_shru1(tt), indc_ave2_shru1(tt)]=max(abs(J_shru1(:)));
    [ii_shru1 jj_shru1]=ind2sub(size(latitude), indc_ave2_shru1(tt));


    %% Correlation map SHRU5
%     nexttile
    subplot(1,2,1)
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

    % find coords of max corr in lat/lon
    [toto, tata]=ind2sub(size(latitude), indc_ave2(tt));    %indc_ave2 is index of max corr
    %%% add location of max point
    plotm(double(latitude(toto, tata)),double(longitude(toto,tata)),'xr','markersize',16,'linewidth',3)

    %     %%% add edges
    %     plotm(latitude_edge_ok, longitude_edge_ok, '.k','markersize',3,'linewidth',1)

    c=colorbar;
    c.Label.String = 'Correlation coefficient';
    caxis(ccc)
    %     title(['2-day averaged SPL - Max corr = ' num2str(R(tt))])
    title([datestr(t_beg_num, 'dd mmm yyyy') ' to ' datestr(t_end_num, 'dd mmm yyyy')]); %, 'fontsize',20,'fontweight', 'bold')

%     ylabel('SHRU5', 'fontsize',30,'fontweight', 'bold')

    dist(tt)=distance(gps_site(1), gps_site(2), double(latitude(toto, tata)),double(longitude(toto, tata)),referenceSphere('Earth'));

    %% Time series SHRU5

    % index within timespan
    SPL_ANL_ok=SPL_ANL_ave2(ind_t_ok_osisaf);
    d_ok=squeeze(d(ind_t_ok_osisaf,ii,jj));
    %    d_ok_interp=interp1(t_osisaf(ind_t_ok_osisaf), d_ok, t(ind_t_ok), 'nearest');
    t_ok=t_osisaf(ind_t_ok_osisaf); %%% common time axis


    ind_no_nan=~isnan(d_ok);


    [SPL_ANL_ok_norm, mu_spl, sigma_spl]=zscore(SPL_ANL_ok(ind_no_nan));
    [d_ok_norm, mu_d, sigma_d]=zscore(d_ok(ind_no_nan));
    d_ok_norm2=(d_ok-mu_d)/sigma_d;

    %%%%%%%%% FINDING CORRELATION VALUES
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

    %     nexttile
    subplot(1,2,2)
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
    legend('ANL (SHRU5)', 'Ice drift')

    % title of the whole figure
    sgtitle(['Ice Drift Correlation for ' num2str(freq_range1) '-' num2str(freq_range2) ' Hz'])
    
    
    
    %% CUTTING VALUES
    %%% plot data over 0.5
    corr_spa_ave2_cut=zeros(119,177,tt);


    % take out anything below 0.45
    for i = 1:119
        for ii = 1:177
            if 1==isnan(corr_spa_ave2(i,ii,tt))
                corr_spa_ave2_cut(i,ii,tt)=nan;
                % do nothing
            elseif (corr_spa_ave2(i,ii,tt))<=0.4
                % change tht index to a nan if less than 0.4
                corr_spa_ave2_cut(i,ii,tt)=nan;
            elseif (corr_spa_ave2(i,ii,tt))>0.4
                corr_spa_ave2_cut(i,ii,tt)=corr_spa_ave2(i,ii,tt);
            end
        end
    end

    % clear so not rewritten every round
    clear r_ind c_ind
    clear dist_map
    clear nans_corr_cut 


    nans_corr_cut = ~isnan(corr_spa_ave2_cut(:,:,tt)); %nan is 0, # is 1
    [r_ind, c_ind,v_corr] = find(nans_corr_cut);    % finds non zeros and saves inds
    % length x whatever lat/long area is shows coverage
    uncut_corr = corr_spa_ave2_cut(r_ind,c_ind,tt);   %

    % find size of boi
    earthellipsoid = referenceSphere('earth','km');
%     area_base = areaquad(double(latitude(1,1)),double(longitude(1,1)),double(latitude(1,2)),double(longitude(1,2)),earthellipsoid,'degrees');
%     area_cov(tt) = areaint(double(latitude(r_ind,c_ind)),double(longitude(r_ind,c_ind)),earthellipsoid)
    
    % count # of pixels active, mulitply by pixel
    %     r_ind*area
    toto=find(corr_spa_ave2_cut(:,:,tt)>0.4) ;
    N_pixel=length(toto);
    pixel_lat=latitude(toto);
    pixel_lon=longitude(toto);
    area_cov(tt) = N_pixel*62.5*62.5    % this is in km!!!!


    % find all the distances
    for iv=1:length(r_ind)
        ri=r_ind(iv);
        ci=c_ind(iv);
        dist_map(iv)=distance(gps_site(1), gps_site(2), double(latitude(ri, ci)),double(longitude(ri,ci)),referenceSphere('Earth'));
    end

    % find index of closest point
    [min_dist(tt), min_ind(tt)]=min(dist_map);
    minlat(tt)=double(latitude(r_ind(min_ind(tt)), c_ind(min_ind(tt))));
    minlon(tt)=double(longitude(r_ind(min_ind(tt)),c_ind(min_ind(tt))));

    % find index of farthest
    [max_dist(tt),max_ind(tt)]=max(dist_map);
    maxlat(tt)=double(latitude(r_ind(max_ind(tt)), c_ind(max_ind(tt))));
    maxlon(tt)=double(longitude(r_ind(max_ind(tt)),c_ind(max_ind(tt))));

    % find index of mid of avg point ( for plotting purposes)
    avg_dist(tt)=mean(dist_map);
    [close_avg_dist(tt), avg_ind(tt)]=min(abs(mean(dist_map)-dist_map));
    avglat(tt)=double(latitude(r_ind(avg_ind(tt)), c_ind(avg_ind(tt))));
    avglon(tt)=double(longitude(r_ind(avg_ind(tt)),c_ind(avg_ind(tt))));

    %%%% make the figure(s) plotting this
    figure
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

    %%% plot CUTOFF data
    surfm(double(latitude), double(longitude), squeeze(corr_spa_ave2_cut(:,:,tt)))

    %%% add mooring
    plotm(gps_site(1),gps_site(2),'xk','markersize',16,'linewidth',3)
    [toto, tata]=ind2sub(size(latitude), indc_ave2(tt));
    %%% add location of max corr point
    plotm(double(latitude(toto, tata)),double(longitude(toto,tata)),'xr','markersize',10,'linewidth',3)

    %%% location of far point
    plotm(maxlat(tt),maxlon(tt),'.m','markersize',20,'linewidth',3)
   
    %%% location of min point
    plotm(minlat(tt),minlon(tt),'.c','markersize',20,'linewidth',3)

    %%% location of avg point
    plotm(avglat(tt),avglon(tt),'dg','markersize',10,'linewidth',3)%'MarkerFaceColor','#D95319')
   
    %     %%% add edges
    %     plotm(latitude_edge_ok, longitude_edge_ok, '.k','markersize',3,'linewidth',1)

    c=colorbar;
    c.Label.String = 'Correlation coefficient';
    caxis([0.4 0.7])
    title(['Correlation over 0.45 for ' num2str(freq_range1) '-' num2str(freq_range2) ' Hz'])
    ylabel([datestr(t_beg_num, 'dd mmm yyyy') ' to ' datestr(t_end_num, 'dd mmm yyyy')])
    xlabel(['MinDist:' num2str(min_dist(tt)) ' MaxDist:' num2str(max_dist(tt)) ' AvgDist:' num2str(avg_dist(tt))])




    % printer currently OFF
    %%%%%%%%%%%% TURN ON AND OFF PRINTING %%%%%%%%%%%%%%%%%%%%%%%%
%         print(gcf,['./new_figs/spatial_corr_result/jasa_plot/' num2str(freq_range1) '_' ... 
%             num2str(freq_range4) '/spatial_corr_' datestr(t_beg_num, 'yyyymmdd') '-' datestr(t_end_num, 'yyyymmdd')]  ...
%             ,'-dpng')

%       

end
%%%% end of the figure generator %%%%
savestring = ['cutoff_' num2str(freq_range1) '_' num2str(freq_range2) '_icecorr']
    save(savestring,'freq_range1','freq_range2','corr_spa_ave2_cut','area_cov','min_dist','minlat','minlon',...
        'max_dist','maxlat','maxlon','avg_dist','avglat','avglon')

end % of freq rounder

%%%%% Functions! %%%%%%%
dist_shrus=distance(gps_site(1), gps_site(2),gps_site_shru1(1), gps_site_shru1(2),referenceSphere('Earth'))/1000;

dist_corr_shru5=[min(dist(3:loop_end)) mean(dist(3:loop_end)) max(dist(3:loop_end))]/1000;

% dist_corr_shru1=[min(dist_shru1(3:loop_end)) mean(dist_shru1(3:loop_end)) max(dist_shru1(3:loop_end))]/1000

r_shru5=[min(R(3:loop_end)) mean(R(3:loop_end)) max(R(3:loop_end))];

% r_shru1=[min(R_shru1(3:loop_end)) mean(R_shru1(3:loop_end)) max(R_shru1(3:loop_end))]
