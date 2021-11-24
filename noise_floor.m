close all
clear all
clc

load ArchivedPSDcomputation/PSD_CANAPE_SHRU5_914_sw_0.0625_lw_420_osw_0.03125_olw_0.mat

t=datenum(timestamp_wavDataFiles,'yyyymmddHHMMSS');
    
%%
figure
imagesc(t, fPSD,10*log10(squeeze(vPSD_pwelch_kinda(:,:,1))).')
datetick2('x')
axis xy

%%
figure
hold on
for cc=1:4
    subplot(2,2,cc)
    toto=squeeze(vPSD_pwelch_kinda(:,:,cc));
    noise=(toto(1,:)+toto(2,:))/2;
   
    semilogx(fPSD, 10*log10(noise))
    grid on
% ylim([25 45])
    xlabel('Frequency (Hz)')
    ylabel('PSD (dB re 1 \muPa^2 / Hz)')

    title(['Channel ' num2str(cc-1)])

end

%%

noise_hydro_dB=zeros(size(fPSD));

noise_hydro_dB(fPSD<100)=53;
noise_hydro_dB(fPSD>100 & fPSD<1000)=40;
noise_hydro_dB(fPSD>1000)=38;

semilogx(fPSD, noise_hydro_dB)

noise_total_dB=10*log10(noise.' +10.^(noise_hydro_dB/10));

figure
semilogx(fPSD, noise_total_dB)


%% Simple computation
10*log10(10^(28/10)+10^(38/10))


% load ArchivedPSDcomputation/PSD_CANAPE_SHRU5_914_sw_1_lw_200_osw_0_olw_0.mat
% 
% t=datenum(timestamp_wavDataFiles,'yyyymmddHHMMSS');
%     
% %%
% figure
% imagesc(t, fPSD,10*log10(squeeze(vPSD_pwelch_kinda(:,:,1))).')
% datetick2('x')
% axis xy
% 
% %%
% figure
% hold on
% for cc=1:4
%     subplot(2,2,cc)
%     toto=squeeze(vPSD_pwelch_kinda(:,:,cc));
%     noise=toto(1,:);
%    
%     semilogx(fPSD, 10*log10(noise))
%     grid on
% % ylim([25 45])
%     xlabel('Frequency (Hz)')
%     ylabel('PSD (dB re 1 \muPa^2 / Hz)')
% 
%     title(['Channel ' num2str(cc-1)])
% 
% end
% 

