close all
clear all
clc

load ArchivedPSDcomputation/PSD_CANAPE_SHRU5_914_beamform_sw_0.0625_lw_420_osw_0.03125_olw_0_theta_0_45_90_135_180.mat
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
Ntheta=length(theta);
f1=[30 250  50  500 1000];
f2=[80 350  500 1000 2000];

ANL=squeeze(vPSD_kinda(:,:,ind_p,:));
LTSA=vPSD_pwelch_kinda;

%%
figure
for tt=1:Ntheta
    
    ax{tt}=subplot(Ntheta,1,tt);
    imagesc(timestamp_num_spectro, fPSD, 10*log10(squeeze(ANL(:,:,tt))).')
    axis xy
    colorbar
    datetick('x')
    % caxis([50 120])
    hold on
    plot(timestamp_num_ssmi,T_ssmi.icefrac*15, 'r')
end
linkaxes([ax{1} ax{2} ax{3} ax{4} ax{5}], 'xy')

figure
for tt=1:Ntheta    
    ax{tt}=subplot(Ntheta,1,tt);
    imagesc(timestamp_num_spectro, fPSD, 10*log10(squeeze(LTSA(:,:,tt))).')
    axis xy
    colorbar
    datetick('x')
    % caxis([50 120])
    hold on
    plot(timestamp_num_ssmi,T_ssmi.icefrac*15, 'r')
end
linkaxes([ax{1} ax{2} ax{3} ax{4} ax{5}], 'xy')


figure
ind=1;
for tt=[1,3]
    
    ax{ind}=subplot(2,1,ind);
    imagesc(timestamp_num_spectro, fPSD, 10*log10(squeeze(ANL(:,:,tt))).')
    axis xy
    colorbar
    datetick('x')
    % caxis([50 120])
    hold on
    plot(timestamp_num_ssmi,T_ssmi.icefrac*15, 'r')
    ind=ind+1;
    ylim([0 500])
end
linkaxes([ax{1} ax{2}], 'xy')

subplot(2,1,1);
title('0 deg (up)')
subplot(2,1,2);
title('90 deg (horizontal)')


%%
Nf=length(f1);
Nt=size(vPSD_pwelch_kinda,1);
SPL_ANL=zeros(Nt,Nf, Ntheta);


for ff=1:Nf
    SPL_ANL(:,ff,:)=10*log10(sum(ANL(:,fPSD>f1(ff) & fPSD<f2(ff),:),2)*(f2(ff)-f1(ff)));
    SPL_raw(:,ff,:)=10*log10(sum(LTSA(:,fPSD>f1(ff) & fPSD<f2(ff),:),2)*(f2(ff)-f1(ff)));
end

for ff=1:Nf
    figure
    subplot(211)
    for tt=1:Ntheta
        plot(1:Nt,squeeze(SPL_ANL(:,ff,tt)))
        hold on
    end
    grid on
    datetick('x')
    title(['Ambient noise level in [' num2str(f1(ff)) ' - ' num2str(f2(ff)) '] Hz'])

    subplot(212)
    for tt=1:Ntheta
        plot(1:Nt,squeeze(SPL_raw(:,ff,tt)))
        hold on
    end
    grid on
    datetick('x')
    title(['Raw sound level in [' num2str(f1(ff)) ' - ' num2str(f2(ff)) '] Hz'])
    legend('0 deg (up)','45 deg', '90 deg (horizontal)','135 deg','180 deg (down)')
end

%%
save ANL_SHRU5_beamform timestamp_num_spectro SPL_ANL SPL_raw  ind_p f1 f2 Nf Nt Ntheta theta ...
    T_ssmi timestamp_num_ssmi T_ecmwf timestamp_num_ecmwf T_smos timestamp_num_smos T_temp timestamp_num_temp 


%%
ff=2;
figure
for tt=1:Ntheta
    plot(timestamp_num_spectro,squeeze(SPL_ANL(:,ff,tt)))
    hold on
end
grid on
datetick2('x')
title(['Ambient noise level in [' num2str(f1(ff)) ' - ' num2str(f2(ff)) '] Hz'])

legend('0 deg (up)','45 deg', '90 deg (horizontal)','135 deg','180 deg (down)')

%
toto_up=squeeze(SPL_ANL(:,2,1));
toto_hor=squeeze(SPL_ANL(:,2,3));

tutu=find(toto_up<toto_hor-10);

toto_up(tutu)=NaN;
toto_hor(tutu)=NaN;


figure
plot(timestamp_num_spectro, toto_up)
hold on
plot(timestamp_num_spectro, toto_hor)
legend('up', 'hor')
datetick2('x')

%%
load ANL_SHRU5_beamform
ff=2;

figure
tt=1;
plot(timestamp_num_spectro,squeeze(SPL_ANL(:,ff,tt)))
hold on
tt=3;
plot(timestamp_num_spectro,squeeze(SPL_ANL(:,ff,tt)))
hold on
grid on
datetick2('x')
title(['Ambient noise level in [' num2str(f1(ff)) ' - ' num2str(f2(ff)) '] Hz'])
legend('0 deg (up)', '90 deg (horizontal)')

figure
tt=1;
plot(timestamp_num_spectro,squeeze(SPL_ANL(:,ff,tt)), 'linewidth', 2)
hold on
tt=3;
plot(timestamp_num_spectro,squeeze(SPL_ANL(:,ff,tt)), 'linewidth', 2)
hold on
grid on
title(['Ambient noise level in [' num2str(f1(ff)) ' - ' num2str(f2(ff)) '] Hz'])
datetick2('x')
legend('0 deg (up)', '90 deg (horizontal)')

diff_anl=squeeze(SPL_ANL(:,ff,3)-SPL_ANL(:,ff,1));
diff_anl2=diff_anl;
diff_anl2(diff_anl2<8)=NaN;
diff_anl2(diff_anl2>12)=NaN;
figure
plot(timestamp_num_spectro,diff_anl)
% hold on
% plot(timestamp_num_spectro,diff_anl2)
grid on
xlim([datenum(2016,11,1) datenum(2017,8,1)])
datetick('x', 'keeplimits')
ylim([-5 15])
title('ANL_{hor}-ANL_{up}')