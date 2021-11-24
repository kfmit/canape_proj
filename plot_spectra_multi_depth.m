close all
clear all
clc

load ArchivedPSDcomputation/PSD_CANAPE_SHRU5_914_sw_0.0625_lw_420_osw_0.03125_olw_0.mat
ind_cut=34;

% load ArchivedPSDcomputation/PSD_sw_1_lw_200_osw_0_olw_0.mat
% ind_cut=1;

% load ArchivedPSDcomputation/PSD_sw_2_lw_300_osw_1_olw_0.mat
% ind_cut=1;

% load ArchivedPSDcomputation/PSD_sw_60_lw_3600_osw_30_olw_0.mat
% ind_cut=1;

vPSD_kinda(1:ind_cut,:,:,:)=[];
vPSD_pwelch_kinda(1:ind_cut,:,:)=[];
timestamp_wavDataFiles(1:ind_cut)=[];

T = readtable('auxData_SHRU5/variables_ASI-SSMI_SHRU5.csv');
timestamp_num_env=datenum(num2str(T.timestamp),'yyyymmddHHMMSS');

timestamp_num_spectro=datenum(timestamp_wavDataFiles,'yyyymmddHHMMSS');


%%
color={'r', 'g', 'k', 'b'};

ind_p=3;

figure
for ind_capt=1:4;

    f1=50;
    f2=500;

    ANL=squeeze(vPSD_kinda(:,:,ind_p,ind_capt));

    LTSA=squeeze(vPSD_pwelch_kinda(:,:,ind_capt));

    SPL_ANL=10*log10(sum(ANL(:,fPSD>f1 & fPSD<f2),2)*(f2-f1));
    SPL_raw=10*log10(sum(LTSA(:,fPSD>f1 & fPSD<f2),2)*(f2-f1));
    %% Define ice / no-ice

    p_no_ice=15;
    ind_no_ice = find(T.icefrac<p_no_ice);
    [~, toto]=max(abs(diff(ind_no_ice)));


    p_ice=85;
    ind_ice=find(T.icefrac>p_ice);


    beg_ice=timestamp_num_env(ind_ice(1));
    end_ice=timestamp_num_env(ind_ice(end));
    ind_acoust_ice=find(timestamp_num_spectro>beg_ice & timestamp_num_spectro<end_ice);


    end_no_ice_acoust=find(timestamp_num_spectro<timestamp_num_env(ind_no_ice(toto-1)));
    end_no_ice_acoust=end_no_ice_acoust(end);
    beg_no_ice_acoust=find(timestamp_num_spectro>timestamp_num_env(ind_no_ice(toto+1)));
    beg_no_ice_acoust=beg_no_ice_acoust(1);
    ind_acoust_no_ice=[1:end_no_ice_acoust beg_no_ice_acoust:length(timestamp_num_spectro)];



    % figure
    % p1=subplot(211);
    % plot(timestamp_num_env,T.icefrac)
    % hold on
    % plot(timestamp_num_env(ind_no_ice),T.icefrac(ind_no_ice), 'o')
    % hold on
    % plot(timestamp_num_env(ind_ice),T.icefrac(ind_ice), '*')
    % datetick('x')
    % 
    % p2=subplot(212);
    % plot(timestamp_num_spectro,SPL_ANL)
    % hold on
    % plot(timestamp_num_spectro(ind_acoust_no_ice),SPL_ANL(ind_acoust_no_ice), 'o')
    % hold on
    % plot(timestamp_num_spectro(ind_acoust_ice),SPL_ANL(ind_acoust_ice), '*')
    % grid on
    % datetick('x')
    % 
    % linkaxes([p1,p2],'x')




    %% Plot spectra

    p_plot=[1 10 25 50 75 90 99];

    perc_ANL_no_ice=prctile(ANL(ind_acoust_no_ice,:),p_plot, 1);
    perc_ANL_ice=prctile(ANL(ind_acoust_ice,:),p_plot, 1);
    perc_LTSA_no_ice=prctile(LTSA(ind_acoust_no_ice,:),p_plot, 1);
    perc_LTSA_ice=prctile(LTSA(ind_acoust_ice,:),p_plot, 1);



    p1=subplot(221);
    hold on
    semilogx(fPSD, 10*log10(perc_ANL_no_ice), 'Color', color{ind_capt})
    grid on
    title('No ice - ANL')
%     legend('1 %', '10 %', '25 %', '50 %', '75 %', '90 %', '99 %')
    xlabel('Frequency (Hz)')
    ylabel('dB / Hz ref 1 \muPa')

    p2=subplot(222);
    hold on
    semilogx(fPSD, 10*log10(perc_ANL_ice), 'Color', color{ind_capt})
    grid on
    title('Ice covered - ANL')
%     legend('1 %', '10 %', '25 %', '50 %', '75 %', '90 %', '99 %')
    xlabel('Frequency (Hz)')
    ylabel('dB / Hz ref 1 \muPa')

    p3=subplot(223);
    hold on
    semilogx(fPSD, 10*log10(perc_LTSA_no_ice), 'Color', color{ind_capt})
    grid on
    title('No ice - LTSA')
%     legend('1 %', '10 %', '25 %', '50 %', '75 %', '90 %', '99 %')
    xlabel('Frequency (Hz)')
    ylabel('dB / Hz ref 1 \muPa')

    p4=subplot(224);
    hold on
    semilogx(fPSD, 10*log10(perc_LTSA_ice), 'Color', color{ind_capt})
    grid on
    title('Ice covered - LTSA')
%     legend('1 %', '10 %', '25 %', '50 %', '75 %', '90 %', '99 %')
    xlabel('Frequency (Hz)')
    ylabel('dB / Hz ref 1 \muPa')

    linkaxes([p1,p2,p3,p4],'xy')
    ylim([30 90])
end
