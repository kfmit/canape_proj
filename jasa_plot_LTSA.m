close all
clear all
clc

load ANL_SHRU5.mat

load ArchivedPSDcomputation/PSD_CANAPE_SHRU5_914_sw_0.0625_lw_420_osw_0.03125_olw_0.mat
ind_cut=34;

vPSD_kinda(1:ind_cut,:,:,:)=[];
vPSD_pwelch_kinda(1:ind_cut,:,:)=[];
timestamp_wavDataFiles(1:ind_cut)=[];

T = readtable('auxData_SHRU5/variables_ASI-SSMI_SHRU5.csv');
timestamp_num_env=datenum(num2str(T.timestamp),'yyyymmddHHMMSS');

T2= readtable('auxData_SHRU5/variables_ECMWF_SHRU5.csv');
timestamp_num_ecmwf=datenum(num2str(T2.timestamp),'yyyymmddHHMMSS');

T3= readtable('auxData_SHRU5/variables_SMOS_SHRU5.csv');
timestamp_num_smos=datenum(num2str(T3.timestamp),'yyyymmddHHMMSS');

T4= readtable('auxData_SHRU5/variables_temp_SHRU5.csv');
timestamp_num_temp=datenum(num2str(T4.timestamp),'yyyymmddHHMMSS');

timestamp_num_spectro=datenum(timestamp_wavDataFiles,'yyyymmddHHMMSS');
%%
ind_p=3;
ind_capt=4;
f1=250;
f2=350;

ANL=squeeze(vPSD_kinda(:,:,ind_p,ind_capt));

LTSA=squeeze(vPSD_pwelch_kinda(:,:,ind_capt));


%%

figure

subplot(311)
imagesc(timestamp_num_spectro, fPSD, 10*log10(LTSA).')
axis xy
% colormap('jet')
% c=colorbar;
% c.Label.String = 'dB re 1\muPa^2';
% datetick('x', 'mmm yyyy')
caxis([20 110])
hold on
plot(timestamp_num_env,T.icefrac*15, 'r', 'linewidth', 2)
title('Power Spectral Density')
ylabel('Frequency (Hz)')
% xlim([timestamp_num_spectro(1) timestamp_num_spectro(end)])
xlim([datenum([2016,10,1]), datenum([2017,08,1])])
xticks([datenum([2016,10,1]) datenum([2016,12,1]) datenum([2017,2,1]) datenum([2017,4,1]) datenum([2017,6,1]) datenum([2017,8,1])])
xticklabels({'Oct 2016','Dec 2016','Feb 2017','Apr 2017','Jun 2017','Aug 2017'})



subplot(312)
imagesc(timestamp_num_spectro, fPSD, 10*log10(ANL).')
axis xy
% colormap('jet')

caxis([20 110])
hold on
plot(timestamp_num_env,T.icefrac*15, 'r', 'linewidth', 2)
title('Ambient Noise Level')
ylabel('Frequency (Hz)')

% xlim([timestamp_num_spectro(1) timestamp_num_spectro(end)])
xlim([datenum([2016,10,1]), datenum([2017,08,1])])
xticks([datenum([2016,10,1]) datenum([2016,12,1]) datenum([2017,2,1]) datenum([2017,4,1]) datenum([2017,6,1]) datenum([2017,8,1])])
xticklabels({'Oct 2016','Dec 2016','Feb 2017','Apr 2017','Jun 2017','Aug 2017'})

c=colorbar('southoutside')
c.Label.String = 'dB re 1\muPa^2 / Hz';


%%

ind_f=5;

subplot(313)
% yyaxis left
plot(timestamp_num_spectro, SPL_raw(:,ind_f))
hold on
plot(timestamp_num_spectro, SPL_ANL(:,ind_f))
ylabel('Sound level (dB re 1\muPa^2 / Hz)')
legend('PSD_{300}', 'ANL_{300}')

% yyaxis right
% hold on
% plot(timestamp_num_env,T.icefrac, 'r', 'linewidth', 2)

grid on
xlim([datenum([2016,10,1]), datenum([2017,08,1])])
xticks([datenum([2016,10,1]) datenum([2016,12,1]) datenum([2017,2,1]) datenum([2017,4,1]) datenum([2017,6,1]) datenum([2017,8,1])])
xticklabels({'Oct 2016','Dec 2016','Feb 2017','Apr 2017','Jun 2017','Aug 2017'})

%%
figure
plot(timestamp_num_spectro, SPL_raw(:,ind_f),'linewidth',2)
hold on
plot(timestamp_num_spectro, SPL_ANL(:,ind_f),'linewidth',2)

% yyaxis right
% hold on
% plot(timestamp_num_env,T.icefrac, 'r', 'linewidth', 2)

grid on
% xlim([datenum([2016,10,1]), datenum([2017,08,1])])
% xticks([datenum([2016,10,1]) datenum([2016,12,1]) datenum([2017,2,1]) datenum([2017,4,1]) datenum([2017,6,1]) datenum([2017,8,1])])
% xticklabels({'Oct 2016','Dec 2016','Feb 2017','Apr 2017','Jun 2017','Aug 2017'})

datetick2('x', 'dd mmm yyyy', 'keeplimits')

