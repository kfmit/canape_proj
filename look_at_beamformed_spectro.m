close all
clear all
clc

load ./ArchivedPSDcomputation/PSD_CANAPE_SHRU5_914_beamform_sw_0.0625_lw_420_osw_0.03125_olw_0_theta_0_45_90_135_180_old.mat
time_num=datenum(timestamp_wavDataFiles,'yyyymmddHHMMSS');

%%%% correct a stupid bug in LTSA_beamform computation
time_num=unique(time_num);


tt=1; %%% angle indice
pp=3; %%% percentile indice

LTSA=10*log10(squeeze(vPSD_pwelch_kinda(:,:,tt))).';
ANL=10*log10(squeeze(vPSD_kinda(:,:,pp,tt))).';

figure
imagesc(time_num, fPSD, LTSA)
axis xy
datetick2('x')

figure
imagesc(time_num, fPSD, ANL)
axis xy
datetick2('x')