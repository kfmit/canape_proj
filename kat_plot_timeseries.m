%% kat's code to plot a couple of timeseries next to each other

close all
clear all
clc

addpath('/home/kfung/Downloads/CANAPE/mat_files/');


load ANL_SHRU5_bigfreqs.mat
load sunrise_sunset_2017_feb_april.mat
ff=5; % this picks the freq

%% originals
figure
subplot(211)
yyaxis left
plot(timestamp_num_spectro,SPL_ANL(:,ff))
% yyaxis right
% plot(timestamp_num_ssmi,T_ssmi.icefrac)
ylim([-10 110])
grid on
datetick('x')
title(['Ambient noise level for [' num2str(f1(ff)) ' - ' num2str(f2(ff)) ']  Hz'])


subplot(212)
yyaxis left
plot(timestamp_num_spectro,SPL_raw(:,ff))
yyaxis right
plot(timestamp_num_ssmi,T_ssmi.icefrac)
ylim([-10 110])
grid on
datetick('x','mmmyy')
title(['Raw data for [' num2str(f1(ff)) ' - ' num2str(f2(ff)) ']  Hz'])
%% new figure

freq_ind = [10 20 30];  %this is for 500 1000 1500

figure
tiledlayout(3,1)

for i=1:3
    ff = freq_ind(i);
    nexttile
%     yyaxis left
    plot(timestamp_num_spectro,SPL_ANL(:,ff))
%     yyaxis right
%     plot(timestamp_num_ssmi,T_ssmi.icefrac)
    ylim([50 110])
    yticks([50  75 100])
    grid on
    datetick('x')
    title(['[' num2str(f1(ff)) ' - ' num2str(f2(ff)) ']  Hz'])
end

sgtitle('Ambient Noise Level Time Series')

%% couple weeks overlay

freq_ind = [10 20 30];  %this is for 500 1000 1500

figure

date_week = [736878:1:736885];
s_date_ind = find(736878 < timestamp_num_spectro);
e_date_ind = find(736885 < timestamp_num_spectro);
date_ind = [s_date_ind(1):1:e_date_ind(1)];

for i=1:3
    ff = freq_ind(i);
%     yyaxis left
    plot(timestamp_num_spectro(date_ind),SPL_ANL(date_ind,ff))
%     yyaxis right
%     plot(timestamp_num_ssmi,T_ssmi.icefrac)
%     title(['[' num2str(f1(ff)) ' - ' num2str(f2(ff)) ']  Hz'])
    hold on
end
hold off
xlim([736878 736885])
% ylim([50 110])
ylabel('SPL ANL')
legend('500 Hz','1000 Hz','1500 Hz','Location','best')

% datetick('x',1)
xticks(date_week)
xticklabels({'02JUL2017', '03JUL2017', '04JUL2017','05JUL2017','06JUL2017','07JUL2017','08JUL2017','09JUL2017',})
xlabel('Date')
grid on
title('Ambient Noise Level Time Series')
legend('500 Hz','1000 Hz','1500 Hz','Location','northwest')
