%% Plot ANL and ice drift
clear all
close all
clc

addpath('/home/kfung/Downloads/CANAPE/mat_files/')
addpath('/home/kfung/Downloads/CANAPE/spatial_mat_files/')
addpath('/home/kfung/Downloads/CANAPE/new_figs/')

freq_array1=[40 250 450 900 1250];
freq_array2=[60 350 550 1100 1750];

figure('visible','on');
tiledlayout(4,2) % plots all
%     tiledlayout(1,2) % plots one

for freqi=2:length(freq_array1)
    freq_range1 = freq_array1(freqi);
    freq_range2 = freq_array2(freqi);
    filename1 = ['spatial_cor_results_interp_new_shru1_' num2str(freq_range1) '_' num2str(freq_range2) '.mat'];
    filename2 = ['spatial_cor_results_interp_new_' num2str(freq_range1) '_' num2str(freq_range2) '.mat'];

    %%%% comes from spatial_corr_analysis_loop_interp_SHRU5
    %% LOAD SHRU5 Data
    % load spatial_cor_results_interp_new.mat and spatial_cor_results_interp_new_1250_1750.mat

    % NEW SHRU5 name
    load(filename2);
    % save the specific variables to new names!
    corr_spa_ave2_shru5_2 = corr_spa_ave2;
    SPL_ANL_ave2_SHRU5_2 = SPL_ANL_ave2;



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


    %% prepare plot

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Loop for generating figure %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURE VISIBILITY HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     figure('visible','on');
    % %     tiledlayout(5,2) % plots all
    %     tiledlayout(1,2) % plots one
    %%%%%%%%%%%%%%%%%%%%%%%%%%5%% FIGURE VISIBILITY HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for tt=3:11, actual loop length is 11
    loop_end = 3; % change to 10 bc fuck it
    for tt=3:loop_end  % loop 1 long
        %     for tt=1
        t_beg_num=date_loop(1,tt);
        t_end_num=date_loop(2,tt);
        ind_t_ok_osisaf=find(t_osisaf >= t_beg_num & t_osisaf <= t_end_num);
        ind_t_raw=find(t >= t_beg_num & t <= t_end_num);

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


        %     figure

        %     h = get(0,'children');
        %     scrsz = get(0,'ScreenSize');
        %     set(h,'Position',[scrsz(1) scrsz(2) scrsz(3)/2 floor(scrsz(4)*0.66)]

        J=squeeze(corr_spa_ave2_shru5_2(:,:,tt));
        J(tutu)=NaN;
        %    figure, imagesc(J)
        [cmax(tt), indc_ave2(tt)]=max(abs(J(:)));
        [ii jj]=ind2sub(size(latitude), indc_ave2(tt));


        %%%%%%%%%%%%%%%%%% SHRU5: BELOW  %%%%%%%%%%%%%%%%%%%

        %% Time series SHRU5, FREQ 1-2
        % ind_t_ok_osiaf gives us the matching date lineups
        % splp_anl_ave is the average per day

        SPL_ANL_ok=SPL_ANL_ave2_SHRU5_2(ind_t_ok_osisaf);
        d_ok=squeeze(d(ind_t_ok_osisaf,ii,jj));
        %    d_ok_interp=interp1(t_osisaf(ind_t_ok_osisaf), d_ok, t(ind_t_ok), 'nearest');
        t_ok=t_osisaf(ind_t_ok_osisaf); %%% common time axis
        % need to find index in spl anl raw that corresponds to start and end of
        % tok then plot(t(ind,eind),spl...

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

        %     subplot(221)
        nexttile
        plot(t_ok-t_ok(1),(SPL_ANL_ok-mu_spl)/sigma_spl, 'Color', [0    0.4470    0.7410], 'linewidth',2)

        hold on
        for bb=1:nb
            plot(t_ok(ind_b{bb})-t_ok(1), block{bb}, '-', 'Color', [0.8500    0.3250    0.0980], 'linewidth',2)
            hold on
        end

        plot(t_ok(ind_no_nan)-t_ok(1), d_ok_norm, 'o', 'Color', [0.8500    0.3250    0.0980], 'linewidth',2)
        title(['R_{max}=' num2str(R(tt))])
        ylabel('Z-score of data')
        xlabel('Days')
        grid on
        ylim([-4 4])

        %%
        %%%%%%%%%%%%%%%%%% plotting the raw ice and ANL %%%%%%%%%%%%%%%%%%%%%%%%
        %     subplot(222)
        nexttile

        %     plot(t_ok-t_ok(1),(SPL_ANL_ok-mu_spl)/sigma_spl, 'Color', [0    0.4470    0.7410], 'linewidth',2)
        yyaxis left
        % plot(t_osisaf,SPL_ANL_ave2_SHRU5_2, 'Color', [0    0.4470    0.7410], 'linewidth',2)
        % plot(t_ok-t_ok(1),SPL_ANL_ok, 'Color', [0    0.4470    0.7410], 'linewidth',2)
        plot(t(ind_t_raw),SPL_ANL(ind_t_raw,freqi), 'Color', [0    0.4470    0.7410], 'linewidth',2)
        ylabel("ANL dB")

        yyaxis right

        % plot(t_ok(ind_no_nan)-t_ok(1), d_ok, '-', 'Color', [0.8500    0.3250    0.0980], 'linewidth',2)
        plot(t_osisaf(ind_t_ok_osisaf),d(ind_t_ok_osisaf,ii,jj), '-', 'Color', [0.8500    0.3250    0.0980], 'linewidth',2)
        title([datestr(t_beg_num, 'dd mmm yyyy') ' to ' datestr(t_end_num, 'dd mmm yyyy')])

        %     title(['R_{max}=' num2str(R(tt))])
        ylabel('Ice drift (km)')
        datetick('x',1)
        grid on

        % title of the whole figure
        %     sgtitle(['Ice Drift and ANL for ' num2str(freq_range1) '-' num2str(freq_range2) ' Hz'])
        sgtitle(['Z-score and Raw Ice Drift and ANL for 300-1500 Hz'])

        % printer currently OFF
        %%%%%%%%%%%% TURN ON AND OFF PRINTING %%%%%%%%%%%%%%%%%%%%%%%%
        %         print(gcf,['./new_figs/spatial_corr_result/jasa_plot/' num2str(freq_range1) '_' ...
        %             num2str(freq_range4) '/spatial_corr_' datestr(t_beg_num, 'yyyymmdd') '-' datestr(t_end_num, 'yyyymmdd')]  ...
        %             ,'-dpng')

    end
    lg=legend(['ANL ' num2str(freq_range1) '-' num2str(freq_range2)], 'Ice drift');
    lg.Location = 'northeastoutside';
end

%% total ice drift
figure
plot(t_osisaf,d(:,ii,jj))
datetick('x')
xlabel('Date')
ylabel('Daily Ice Drift (km)')
title(['Ice Drift from ' datestr(t_osisaf(1), 'dd mmm yyyy') ' to ' datestr(t_osisaf(end), 'dd mmm yyyy')])

%%%% end of the figure generator %%%%


%%%%% Functions! %%%%%%%
% dist_shrus=distance(gps_site(1), gps_site(2),gps_site_shru1(1), gps_site_shru1(2),referenceSphere('Earth'))/1000
%
% dist_corr_shru5=[min(dist(3:loop_end)) mean(dist(3:loop_end)) max(dist(3:loop_end))]/1000
%
% dist_corr_shru1=[min(dist_shru1(3:loop_end)) mean(dist_shru1(3:loop_end)) max(dist_shru1(3:loop_end))]/1000
%
% r_shru5=[min(R(3:loop_end)) mean(R(3:loop_end)) max(R(3:loop_end))]
%
% r_shru1=[min(R_shru1(3:loop_end)) mean(R_shru1(3:loop_end)) max(R_shru1(3:loop_end))]
