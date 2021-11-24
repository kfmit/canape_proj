close all
clear all
clc

plot_env=1;
plot_acoust=1;

color_shru=['r', 'g', ' ','b', 'k'];

if plot_env
%% SHRU 1
T{1} = readtable('auxData_SHRU1/variables_ASI-SSMI_SHRU1.csv');
timestamp_num_env{1}=datenum(num2str(T{1}.timestamp),'yyyymmddHHMMSS');

T2{1}= readtable('auxData_SHRU1/variables_ECMWF_SHRU1.csv');
timestamp_num_ecmwf{1}=datenum(num2str(T2{1}.timestamp),'yyyymmddHHMMSS');

T3{1}= readtable('auxData_SHRU1/variables_SMOS_SHRU1.csv');
timestamp_num_smos{1}=datenum(num2str(T3{1}.timestamp),'yyyymmddHHMMSS');

T4{1}= readtable('auxData_SHRU1/variables_temp_SHRU1.csv');
timestamp_num_temp{1}=datenum(num2str(T4{1}.timestamp),'yyyymmddHHMMSS');

%% SHRU 2
T{2} = readtable('auxData_SHRU2/variables_ASI-SSMI_SHRU2.csv');
timestamp_num_env{2}=datenum(num2str(T{2}.timestamp),'yyyymmddHHMMSS');

T2{2}= readtable('auxData_SHRU2/variables_ECMWF_SHRU2.csv');
timestamp_num_ecmwf{2}=datenum(num2str(T2{2}.timestamp),'yyyymmddHHMMSS');

T3{2}= readtable('auxData_SHRU2/variables_SMOS_SHRU2.csv');
timestamp_num_smos{2}=datenum(num2str(T3{2}.timestamp),'yyyymmddHHMMSS');

T4{2}= readtable('auxData_SHRU2/variables_temp_SHRU2.csv');
timestamp_num_temp{2}=datenum(num2str(T4{2}.timestamp),'yyyymmddHHMMSS');



%% SHRU 4
T{4} = readtable('auxData_SHRU4/variables_ASI-SSMI_SHRU4.csv');
timestamp_num_env{4}=datenum(num2str(T{4}.timestamp),'yyyymmddHHMMSS');

T2{4}= readtable('auxData_SHRU4/variables_ECMWF_SHRU4.csv');
timestamp_num_ecmwf{4}=datenum(num2str(T2{4}.timestamp),'yyyymmddHHMMSS');

T3{4}= readtable('auxData_SHRU4/variables_SMOS_SHRU4.csv');
timestamp_num_smos{4}=datenum(num2str(T3{4}.timestamp),'yyyymmddHHMMSS');

T4{4}= readtable('auxData_SHRU4/variables_temp_SHRU4.csv');
timestamp_num_temp{4}=datenum(num2str(T4{4}.timestamp),'yyyymmddHHMMSS');

%% SHRU 5

T{5} = readtable('auxData_SHRU5/variables_ASI-SSMI_SHRU5.csv');
timestamp_num_env{5}=datenum(num2str(T{5}.timestamp),'yyyymmddHHMMSS');

T2{5}= readtable('auxData_SHRU5/variables_ECMWF_SHRU5.csv');
timestamp_num_ecmwf{5}=datenum(num2str(T2{5}.timestamp),'yyyymmddHHMMSS');

T3{5}= readtable('auxData_SHRU5/variables_SMOS_SHRU5.csv');
timestamp_num_smos{5}=datenum(num2str(T3{5}.timestamp),'yyyymmddHHMMSS');

T4{5}= readtable('auxData_SHRU5/variables_temp_SHRU5.csv');
timestamp_num_temp{5}=datenum(num2str(T4{5}.timestamp),'yyyymmddHHMMSS');


%% Env plots without temp
figure

for ii=1:5
    
    if ii~=3
    
    p1=subplot(411);
    plot(timestamp_num_env{ii},T{ii}.icefrac,color_shru(ii))
    ylabel('Ice Fraction (%)')
    grid on
    hold on
    datetick2('x')
    

    p2=subplot(412);
    plot(timestamp_num_ecmwf{ii}, T2{ii}.W10,color_shru(ii))
    ylabel('Wind speed (m/s)')
    grid on
    hold on
    datetick2('x')

    p3=subplot(413);
    plot(timestamp_num_ecmwf{ii}, T2{ii}.tp,color_shru(ii))
    ylabel('Total precipitation (m) ')
    grid on
    hold on
    datetick2('x')

    p4=subplot(414);
    plot(timestamp_num_smos{ii}, T3{ii}.sea_ice_thickness,color_shru(ii))
    ylabel('Sea ice thickness (m) ')
    grid on
    hold on
    datetick2('x')

    % p6=subplot(616);
    % plot(timestamp_num_temp, T4.T2085, timestamp_num_temp, T4.T2073, timestamp_num_temp, T4.T2050, timestamp_num_temp, T4.sbe3125)
    % ylabel('Underwater temperature (degree C) ')
    % datetick2('x')
    % grid on
    end
end
linkaxes([p1,p2,p3,p4],'x')
xlim([datenum([2016,11,1]) datenum([2017,09,1])])


%% Temp plots


figure
p1=subplot(511);
plot(timestamp_num_temp{1}, T4{1}.T2085, color_shru(1))
hold on
plot(timestamp_num_temp{2}, T4{2}.T2045, color_shru(2))
%%%% no data for shru 4
%%%% no data for shru 5
grid on
title('Depth = 66 m')
datetick2('x')

p2=subplot(512);
plot(timestamp_num_temp{1}, T4{1}.sbe3125, color_shru(1))
hold on
plot(timestamp_num_temp{2}, T4{2}.sbe3123, color_shru(2))
plot(timestamp_num_temp{4}, T4{4}.sbe3078, color_shru(4))
%%%% no data for shru 5
grid on
title('Depth = 141.5 m')
datetick2('x')

p3=subplot(513);
plot(timestamp_num_temp{1}, T4{1}.T2073, color_shru(1))
hold on
plot(timestamp_num_temp{2}, T4{2}.T2074, color_shru(2))
plot(timestamp_num_temp{4}, T4{4}.T2046, color_shru(4))
plot(timestamp_num_temp{5}, T4{5}.T2062, color_shru(5))
grid on
title('Depth = 158 m')
datetick2('x')

p4=subplot(514);
plot(timestamp_num_temp{1}, T4{1}.T2050, color_shru(1))
hold on
%%%% no date for shru 2
plot(timestamp_num_temp{4}, T4{4}.T2049, color_shru(4))
plot(timestamp_num_temp{5}, T4{5}.T2061, color_shru(5))
grid on
title('Depth = 296 m')
datetick2('x')

p5=subplot(515);
plot(timestamp_num_temp{5}, T4{5}.T2075, color_shru(5))
grid on
title('Depth = 446 m')
datetick2('x')


linkaxes([p1,p2,p3,p4,p5],'xy')
xlim([datenum([2016,11,1]) datenum([2017,09,1])])
% xlim([datenum([2017,1,1]) datenum([2017,01,31])])
ylim([-2 2])

end

%% Acoustics

if plot_acoust

ind_p=3;
ind_capt=4;
f1=[50   500   1000];
f2=[500 1000   2000];
Nf=length(f1);
    
%% SHRU 1
ii=1;

load ArchivedPSDcomputation/PSD_CANAPE_SHRU1_903_sw_0.0625_lw_420_osw_0.03125_olw_0.mat
ind_cut=34;

vPSD_kinda(1:ind_cut,:,:,:)=[];
vPSD_pwelch_kinda(1:ind_cut,:,:)=[];
timestamp_wavDataFiles(1:ind_cut)=[];
timestamp_num_spectro{ii}=datenum(timestamp_wavDataFiles,'yyyymmddHHMMSS');
    
ANL=squeeze(vPSD_kinda(:,:,ind_p,ind_capt));
LTSA=squeeze(vPSD_pwelch_kinda(:,:,ind_capt));

for ff=1:Nf
    SPL_ANL{ii}(ff,:)=10*log10(sum(ANL(:,fPSD>f1(ff) & fPSD<f2(ff)),2)*(f2(ff)-f1(ff)));
    SPL_raw{ii}(ff,:)=10*log10(sum(LTSA(:,fPSD>f1(ff) & fPSD<f2(ff)),2)*(f2(ff)-f1(ff)));   
end


%% SHRU 2
ii=2;
load ArchivedPSDcomputation/PSD_CANAPE_SHRU2_906_sw_0.0625_lw_420_osw_0.03125_olw_0.mat
ind_cut=34;

vPSD_kinda(1:ind_cut,:,:,:)=[];
vPSD_pwelch_kinda(1:ind_cut,:,:)=[];
timestamp_wavDataFiles(1:ind_cut)=[];
timestamp_num_spectro{ii}=datenum(timestamp_wavDataFiles,'yyyymmddHHMMSS');
    
ANL=squeeze(vPSD_kinda(:,:,ind_p,ind_capt));
LTSA=squeeze(vPSD_pwelch_kinda(:,:,ind_capt));

for ff=1:Nf
    SPL_ANL{ii}(ff,:)=10*log10(sum(ANL(:,fPSD>f1(ff) & fPSD<f2(ff)),2)*(f2(ff)-f1(ff)));
    SPL_raw{ii}(ff,:)=10*log10(sum(LTSA(:,fPSD>f1(ff) & fPSD<f2(ff)),2)*(f2(ff)-f1(ff)));   
end

%% SHRU 4
ii=4;
load ArchivedPSDcomputation/PSD_CANAPE_SHRU4_905_sw_0.0625_lw_420_osw_0.03125_olw_0.mat
ind_cut=110;

vPSD_kinda(1:ind_cut,:,:,:)=[];
vPSD_pwelch_kinda(1:ind_cut,:,:)=[];
timestamp_wavDataFiles(1:ind_cut)=[];
timestamp_num_spectro{ii}=datenum(timestamp_wavDataFiles,'yyyymmddHHMMSS');
    
ANL=squeeze(vPSD_kinda(:,:,ind_p,ind_capt));
LTSA=squeeze(vPSD_pwelch_kinda(:,:,ind_capt));

for ff=1:Nf
    SPL_ANL{ii}(ff,:)=10*log10(sum(ANL(:,fPSD>f1(ff) & fPSD<f2(ff)),2)*(f2(ff)-f1(ff)));
    SPL_raw{ii}(ff,:)=10*log10(sum(LTSA(:,fPSD>f1(ff) & fPSD<f2(ff)),2)*(f2(ff)-f1(ff)));   
end


%% SHRU 5
ii=5;
load ArchivedPSDcomputation/PSD_CANAPE_SHRU5_914_sw_0.0625_lw_420_osw_0.03125_olw_0.mat
ind_cut=34;

vPSD_kinda(1:ind_cut,:,:,:)=[];
vPSD_pwelch_kinda(1:ind_cut,:,:)=[];
timestamp_wavDataFiles(1:ind_cut)=[];
timestamp_num_spectro{ii}=datenum(timestamp_wavDataFiles,'yyyymmddHHMMSS');
    
ANL=squeeze(vPSD_kinda(:,:,ind_p,ind_capt));
LTSA=squeeze(vPSD_pwelch_kinda(:,:,ind_capt));

for ff=1:Nf
    SPL_ANL{ii}(ff,:)=10*log10(sum(ANL(:,fPSD>f1(ff) & fPSD<f2(ff)),2)*(f2(ff)-f1(ff)));
    SPL_raw{ii}(ff,:)=10*log10(sum(LTSA(:,fPSD>f1(ff) & fPSD<f2(ff)),2)*(f2(ff)-f1(ff)));   
end


%%
figure
ax=[];
for ff=1:Nf
    for ii=1:5
        if ii~=3
            ppp=subplot(Nf,1,ff);
            plot(timestamp_num_spectro{ii},SPL_ANL{ii}(ff,:), color_shru(ii))
            hold on
        end
        ax=[ax, ppp];
    end
    grid on
    hold on
    title(['Ambient noise power in [' num2str(f1(ff)) ' - ' num2str(f2(ff)) '] Hz'])
    datetick2('x')
end
linkaxes(ax,'xy')
xlim([datenum([2016,11,1]) datenum([2017,09,1])])
% xlim([datenum([2017,1,1]) datenum([2017,01,31])])
clear p

figure
for ff=1:Nf
    for ii=1:5
        if ii~=3
            ppp=subplot(Nf,1,ff);
            plot(timestamp_num_spectro{ii},SPL_raw{ii}(ff,:), color_shru(ii))
            hold on
        end
        ax=[ax, ppp];
    end
    grid on
    hold on
    title(['Recorded sound power in [' num2str(f1(ff)) ' - ' num2str(f2(ff)) '] Hz'])
    datetick2('x')
end
linkaxes(ax,'xy')
xlim([datenum([2016,11,1]) datenum([2017,09,1])])

end