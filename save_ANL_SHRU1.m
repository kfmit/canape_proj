close all
clear all
clc

% load ArchivedPSDcomputation/PSD_CANAPE_SHRU5_914_sw_0.0625_lw_420_osw_0.03125_olw_0.mat
load ArchivedPSDcomputation/PSD_CANAPE_SHRU1_903_sw_0.0625_lw_420_osw_0.03125_olw_0.mat
ind_cut=34;

vPSD_kinda(1:ind_cut,:,:,:)=[];
vPSD_pwelch_kinda(1:ind_cut,:,:)=[];
timestamp_wavDataFiles(1:ind_cut)=[];

T_ssmi = readtable('auxData_SHRU5/variables_ASI-SSMI_SHRU5.csv');
timestamp_num_ssmi=datenum(num2str(T_ssmi.timestamp),'yyyymmddHHMMSS');

T_ecmwf= readtable('auxData_SHRU5/variables_ECMWF_SHRU5.csv');
timestamp_num_ecmwf=datenum(num2str(T_ecmwf.timestamp),'yyyymmddHHMMSS');

T_smos= readtable('auxData_SHRU5/variables_SMOS_SHRU5.csv');
timestamp_num_smos=datenum(num2str(T_smos.timestamp),'yyyymmddHHMMSS');

T_temp= readtable('auxData_SHRU5/variables_temp_SHRU5.csv');
timestamp_num_temp=datenum(num2str(T_temp.timestamp),'yyyymmddHHMMSS');

timestamp_num_spectro=datenum(timestamp_wavDataFiles,'yyyymmddHHMMSS');
%%
ind_p=3;
ind_capt=4;
% f1=[30  50  500 1000 250];
% f2=[80 500 1000 2000 350];

%%% NEW Frequencies that I do
f1=[40 450 900 1250 250];
f2=[60 550 1100 1750 350];

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
plot(timestamp_num_ssmi,T_ssmi.icefrac*15, 'r')
title('Ambient noise level')

subplot(212)
imagesc(timestamp_num_spectro, fPSD, 10*log10(LTSA).')
axis xy
colorbar
datetick('x')
% caxis([50 120])
hold on
plot(timestamp_num_ssmi,T_ssmi.icefrac*15, 'r')
title('Raw data')


%%
Nf=length(f1);
Nt=size(vPSD_pwelch_kinda,1);
SPL_ANL=zeros(Nt,Nf);

for ff=1:Nf
    SPL_ANL(:,ff)=10*log10(sum(ANL(:,fPSD>f1(ff) & fPSD<f2(ff)),2)*(f2(ff)-f1(ff)));
    SPL_raw(:,ff)=10*log10(sum(LTSA(:,fPSD>f1(ff) & fPSD<f2(ff)),2)*(f2(ff)-f1(ff)));

    figure
    subplot(211)
    yyaxis left
    plot(timestamp_num_spectro,SPL_ANL(:,ff))
    yyaxis right
    plot(timestamp_num_ssmi,T_ssmi.icefrac)
    ylim([-10 110])
    grid on
    datetick('x')
    title('Ambient noise level')

    subplot(212)
    yyaxis left
    plot(timestamp_num_spectro,SPL_raw(:,ff))
    yyaxis right
    plot(timestamp_num_ssmi,T_ssmi.icefrac)
    ylim([-10 110])
    grid on
    datetick('x','mmmyy')
    title('Raw data')
end

save ANL_SHRU1 timestamp_num_spectro SPL_ANL SPL_raw  ind_p ind_capt f1 f2 Nf Nt ...
    T_ssmi timestamp_num_ssmi T_ecmwf timestamp_num_ecmwf T_smos timestamp_num_smos T_temp timestamp_num_temp 