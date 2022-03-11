%%%% A code to let me put all the spatial corr into one button


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

for i=1:5 %(length(freq_array1)-1)
    freq_range1 = freq_array1(i);
    freq_range2 = freq_array2(i);
    filename1 = ['spatial_cor_results_interp_new_shru1_' num2str(freq_range1) '_' num2str(freq_range2) '.mat'];
    filename2 = ['spatial_cor_results_interp_new_' num2str(freq_range1) '_' num2str(freq_range2) '.mat'];

    load(filename1);
    % save the specific variables to new names!
    corr_spa_ave2_shru1 = corr_spa_ave2;
    SPL_ANL_ave2_SHRU1 = SPL_ANL_ave2;
    gps_site = [72+54.4123/60 , -(159+1.0840/60)];

    %%%% comes from spatial_corr_analysis_loop_interp_SHRU5
    %% LOAD SHRU5 Data
    % load spatial_cor_results_interp_new.mat and spatial_cor_results_interp_new_1250_1750.mat

    % NEW SHRU5 name
    load(filename2);
    % save the specific variables to new names!
    corr_spa_ave2_shru5_2 = corr_spa_ave2;
    SPL_ANL_ave2_SHRU5_2 = SPL_ANL_ave2;

    % load(filename4);
    % % save the specific variables to new names!
    % corr_spa_ave2_shru5_4 = corr_spa_ave2;
    % SPL_ANL_ave2_SHRU5_4 = SPL_ANL_ave2;


    % GPS data
    gps_site = [72+54.4580/60 , -(157+29.2442/60)];
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

    Nloop=size(date_loop,2);

    %% prepare plot map
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
    % for tt=3:11, actual loop length is 11, actual length is 15, but summer
    % months irrelevant
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURE VISIBILITY HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure('visible','on');
    %%%%%%%%%%%%%%%%%%%%%%%%%%5%% FIGURE VISIBILITY HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    loop_end = 11; % change back to 11 after testing
    for tt=3:loop_end  % loop 1 long

        t_beg_num=date_loop(1,tt);
        t_end_num=date_loop(2,tt);
        ind_t_ok_osisaf=find(t_osisaf >= t_beg_num & t_osisaf <= t_end_num);

        %% Read ice edge and type

        t_mid=(t_beg_num+t_end_num)/2;
        ddd_mid=(datestr(t_mid, 'yyyymmdd'));


        h = get(0,'children');
        scrsz = get(0,'ScreenSize');
        set(h,'Position',[scrsz(1) scrsz(2) scrsz(3)/2 floor(scrsz(4)*0.66)])

        J=squeeze(corr_spa_ave2_shru5_2(:,:,tt));
        J(tutu)=NaN;

        %    figure, imagesc(J)
        [cmax(tt), indc_ave2(tt)]=max(abs(J(:)));
        [ii jj]=ind2sub(size(latitude), indc_ave2(tt));


        J_shru1=squeeze(corr_spa_ave2_shru5_2(:,:,tt));
        J_shru1(tutu)=NaN;

        %    figure, imagesc(J)

        % index of max corr 
        [toto, tata]=ind2sub(size(latitude), indc_ave2(tt));


        %%%%%%%%%%%%%%%%%% SHRU5: BELOW  %%%%%%%%%%%%%%%%%%%
        %
        
        % index of max corr 
        [toto, tata]=ind2sub(size(latitude), indc_ave2(tt));

        % save the coords of location
        % useful saveable vars: maxcorr_lat, maxcorr_lon, cmax, dist
        maxcorr_lat(tt)=double(latitude(toto, tata));
        maxcorr_lon(tt)=double(longitude(toto,tata));
        dist(tt)=distance(gps_site(1), gps_site(2), double(latitude(toto, tata)),double(longitude(toto, tata)),referenceSphere('Earth'));        
        %% Time series SHRU5, FREQ 1-2

        SPL_ANL_ok=SPL_ANL_ave2(ind_t_ok_osisaf);
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

    end        % end of loop looking through the months
    % pause inside fig gen to check points have been added

%     % Correlation map SHRU5, FREQ 1-2
%         hold on
%         maph=axesm('MapProjection','lambertstd','MapLatLimit',latlimit,'MapLonLimit',lonlimit);
% 
%         %%% add grid
%         setm(maph,'Grid','on','Glinewidth',1,'PLineLocation',parallel, 'MLineLocation',MLabelLocation);
% 
%         %%% add grid labeling
%         setm(maph,'Fontangle','normal',...
%             'FontSize',12,'fontweight','b',...
%             'MeridianLabel','on',...
%             'MLabelLocation',MLabelLocation,...
%             'MLabelParallel',Mpos,...
%             'ParallelLabel','on',...
%             'PLabelLocation',parallel,...
%             'PLabelMeridian',PLabelMeridian);
% 
%         %%% add land
%         geoshow(maph, land, 'FaceColor',[0.80 0.80 0.80],'EdgeColor',0.30*[1 1 1]);
% 
%         %%% plot data
%         %     surfm(double(latitude), double(longitude), squeeze(corr_spa_ave2(:,:,tt)))
% 
%         %%% add mooring SHRU
%         plotm(gps_site(1),gps_site(2),'xk','markersize',16,'linewidth',3)
% 
%         % IMPORTANT: Location of
%         plotm(maxcorr_lat,maxcorr_lon,'--xr','markersize',16,'linewidth',3)
% 
% %             %%% add edges
% %             plotm(latitude_edge_ok, longitude_edge_ok, '.k','markersize',3,'linewidth',1)
% 
%         c=colorbar;
%         c.Label.String = 'Correlation coefficient';
%         caxis(ccc)
%         %     title(['2-day averaged SPL - Max corr = ' num2str(R(tt))])
%         title([datestr(t_beg_num, 'dd mmm yyyy') ' to ' datestr(t_end_num, 'dd mmm yyyy')], 'fontsize',20,'fontweight', 'bold')
% 
%         ylabel([num2str(freq_range1) '-' num2str(freq_range2)], 'fontsize',30,'fontweight', 'bold')

%         dist(tt)=distance(gps_site(1), gps_site(2), double(latitude(toto, tata)),double(longitude(toto, tata)),referenceSphere('Earth'));

    % printer currently OFF
    %%%%%%%%%%%% TURN ON AND OFF PRINTING %%%%%%%%%%%%%%%%%%%%%%%%
    %         print(gcf,['./new_figs/spatial_corr_result/jasa_plot/' num2str(freq_range1) '_' ...
    %             num2str(freq_range4) '/spatial_corr_' datestr(t_beg_num, 'yyyymmdd') '-' datestr(t_end_num, 'yyyymmdd')]  ...
    %             ,'-dpng')

    savestring = [num2str(freq_range1) '_' num2str(freq_range2) 'corr']
    save(savestring,'freq_range1','freq_range2','maxcorr_lat','maxcorr_lon','cmax','dist')
end     % end of loop going through freqs

%%%% end of the figure generator %%%%


%%%%% Functions! %%%%%%%
dist_shrus=distance(gps_site(1), gps_site(2),gps_site(1), gps_site(2),referenceSphere('Earth'))/1000

dist_corr_shru5=[min(dist(3:loop_end)) mean(dist(3:loop_end)) max(dist(3:loop_end))]/1000

% dist_corr_shru1=[min(dist_shru1(3:loop_end)) mean(dist_shru1(3:loop_end)) max(dist_shru1(3:loop_end))]/1000

r_shru5=[min(R(3:loop_end)) mean(R(3:loop_end)) max(R(3:loop_end))]

% r_shru1=[min(R_shru1(3:loop_end)) mean(R_shru1(3:loop_end)) max(R_shru1(3:loop_end))]
