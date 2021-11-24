close all
clear all
clc





load ArchivedPSDcomputation/PSD_CANAPE_SHRU5_914_sw_1_lw_200_osw_0_olw_0.mat
ind_cut=1; 
ind_capt=4;



vPSD_kinda(1:ind_cut,:,:,:)=[];
vPSD_pwelch_kinda(1:ind_cut,:,:)=[];
timestamp_wavDataFiles(1:ind_cut)=[];


timestamp_num_spectro=datenum(timestamp_wavDataFiles,'yyyymmddHHMMSS');

ind_p=3;
% LTSA=squeeze(vPSD_kinda(:,:,ind_p,ind_capt));
LTSA=squeeze(vPSD_pwelch_kinda(:,:,ind_capt));

%%
clc

N_month=9;
m0=10;

mean_dsp=zeros(N_month,length(fPSD));

for mm=1:N_month
    t0=datenum([2016,m0+mm,1]);
    t1=datenum([2016,m0+mm,31]);
    disp(['from ', datestr(t0), ' to ' datestr(t1)])
    
    ind_t=find(timestamp_num_spectro>=t0 & timestamp_num_spectro<=t1);
    
    LTSA_ok=LTSA(ind_t,:);
    
    mean_dsp(mm,:)=median(LTSA_ok,1);
    
end

%% Figure
co = [         0         0    1.0000 ;
         0    0.5000         0;
    1.0000         0         0;
         0    0.7500    0.7500;
    0.7500         0    0.7500;
    0.7500    0.7500         0;
    0.2500    0.2500    0.2500;
    0.4660    0.6740    0.1880;
    0.8500    0.3250    0.0980];
    
    


figure
for mm=1:7    
    semilogx(fPSD, 10*log10(mean_dsp(mm,:)), 'linewidth',2, 'color', co(mm+2,:))
    hold on
end
% for mm=8:9   
%     semilogx(fPSD, 10*log10(mean_dsp(mm,:)),'--', 'linewidth',2)
%     hold on
% end
grid on
xlim([20 350])
ylim([45 100])
xticks([10 20 30 40 50 100 200 250 350])
xlabel('Frequency (Hz)')
ylabel('Sound Spectrum level (dB re 1 \muPa^2/Hz)')
% legend('November', 'December', 'January', 'February', 'March', 'April', 'May', 'June', 'July')
legend('November', 'December', 'January', 'February', 'March', 'April', 'May')


