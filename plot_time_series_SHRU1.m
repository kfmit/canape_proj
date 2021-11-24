close all
clear all
clc

% load ArchivedPSDcomputation/PSD_CANAPE_SHRU5_914_sw_0.0625_lw_420_osw_0.03125_olw_0.mat
load ArchivedPSDcomputation/PSD_CANAPE_SHRU1_903_sw_0.0625_lw_420_osw_0.03125_olw_0.mat
ind_cut=34;

vPSD_kinda(1:ind_cut,:,:,:)=[];
vPSD_pwelch_kinda(1:ind_cut,:,:)=[];
timestamp_wavDataFiles(1:ind_cut)=[];

T = readtable('auxData_SHRU1/variables_ASI-SSMI_SHRU1.csv');
timestamp_num_env=datenum(num2str(T.timestamp),'yyyymmddHHMMSS');

T2= readtable('auxData_SHRU1/variables_ECMWF_SHRU1.csv');
timestamp_num_ecmwf=datenum(num2str(T2.timestamp),'yyyymmddHHMMSS');

T3= readtable('auxData_SHRU1/variables_SMOS_SHRU1.csv');
timestamp_num_smos=datenum(num2str(T3.timestamp),'yyyymmddHHMMSS');

T4= readtable('auxData_SHRU1/variables_temp_SHRU1.csv');
timestamp_num_temp=datenum(num2str(T4.timestamp),'yyyymmddHHMMSS');

timestamp_num_spectro=datenum(timestamp_wavDataFiles,'yyyymmddHHMMSS');
%%
ind_p=3;
ind_capt=4;
f1=50;
f2=500;

ANL=squeeze(vPSD_kinda(:,:,ind_p,ind_capt));
LTSA=squeeze(vPSD_pwelch_kinda(:,:,ind_capt));

%%
figure
subplot(211)
imagesc(timestamp_num_spectro, fPSD, 10*log10(ANL).')
axis xy
colorbar
datetick('x')
% caxis([50 120])
hold on
plot(timestamp_num_env,T.icefrac*15, 'r')
title('Ambient noise level')

subplot(212)
imagesc(timestamp_num_spectro, fPSD, 10*log10(LTSA).')
axis xy
colorbar
datetick('x')
% caxis([50 120])
hold on
plot(timestamp_num_env,T.icefrac*15, 'r')
title('Raw data')


%%
% SPL=10*log10(mean(ANL(:,fPSD>f1 & fPSD<f2),2));
SPL_ANL=10*log10(sum(ANL(:,fPSD>f1 & fPSD<f2),2)*(f2-f1));
SPL_raw=10*log10(sum(LTSA(:,fPSD>f1 & fPSD<f2),2)*(f2-f1));

figure
subplot(211)
yyaxis left
plot(timestamp_num_spectro,SPL_ANL)
yyaxis right
plot(timestamp_num_env,T.icefrac)
ylim([-10 110])
grid on
datetick('x')
title('Ambient noise level')

subplot(212)
yyaxis left
plot(timestamp_num_spectro,SPL_raw)
yyaxis right
plot(timestamp_num_env,T.icefrac)
ylim([-10 110])
grid on
datetick('x','mmmyy')
title('Raw data')


%%

figure

p1=subplot(611);
plot(timestamp_num_spectro,SPL_ANL)
hold on
plot(timestamp_num_spectro,SPL_raw)
datetick('x')
grid on
ylabel('ANL and SPL (dB re 1 \uPa)')

p2=subplot(612);
plot(timestamp_num_env,T.icefrac)
ylabel('Ice Fraction (%)')
datetick('x')
grid on

p3=subplot(613);
plot(timestamp_num_ecmwf, T2.W10)
ylabel('Wind speed (m/s)')
datetick('x')
grid on

p4=subplot(614);
plot(timestamp_num_ecmwf, T2.tp)
ylabel('Total precipitation (m) ')
datetick('x')
grid on

p5=subplot(615);
yyaxis left
plot(timestamp_num_smos, T3.sea_ice_thickness)
ylabel('Sea ice thickness (m) ')
yyaxis right
plot(timestamp_num_smos, T3.Tsurf)
ylabel('Surface Temp (degree C) ')
datetick('x')
grid on

p6=subplot(616);
plot(timestamp_num_temp, T4.T2085, timestamp_num_temp, T4.T2073, timestamp_num_temp, T4.T2050, timestamp_num_temp, T4.sbe3125)
ylabel('Underwater temperature (degree C) ')
datetick('x')
grid on

linkaxes([p1,p2,p3,p4,p5,p6],'x')
xlim([timestamp_num_spectro(1) timestamp_num_spectro(end)])