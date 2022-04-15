%% kraken trend analysis
% there is a built in function but its kinda weird

clear all
close all
clc

addpath('/home/kfung/Downloads/CANAPE/kraken_data/singleduct_900km_5SD_attn/');
addpath('/home/kfung/Downloads/CANAPE/kraken_data/5mSD_900km_166RD_attn/');
addpath('/home/kfung/Downloads/CANAPE/kraken_data/doubleduct_900km_5SD_ attn/');
% addpath('/home/kfung/Downloads/CANAPE/kraken_data/constant_tester/');
names_dir = dir('5mSD_900km_166RD_attn/');
% names_dir = dir('doubleduct_900km_5SD_ attn/');

actual_order=[7 5 6 3 4];
figure(1)
k = 1000;
for i=1:length(actual_order)
    act_ind = actual_order(i)
    disp(names_dir(act_ind).name)
fig=openfig(names_dir(act_ind).name,'invisible');  %'invisible'
axObjs = fig.Children;
dataObjs = axObjs.Children;
x = dataObjs(1).XData;
y  = dataObjs(1).YData;
R = movmean(x,k);
TL = movmean(y,k);
figure(1)
plot(R,TL)
hold on
end

% for i=3:length(names_dir)-1
%     disp(names_dir(i))
% fig=open(names_dir(i).name);
% axObjs = fig.Children;
% dataObjs = axObjs.Children;
% x(:,i) = dataObjs(1).XData;
% y(:,i)  = dataObjs(1).YData;
% end
% 
% %%
% figure
% k = 1000;
% actual_order=[7 5 6 3 4];
% for ii=1:length(actual_order)
% ind =actual_order(ii);
% R = movmean(x(:,ind),k);
% TL = movmean(y(:,ind),k);
% plot(R,TL)
% hold on
% end

set(gca,'YDir','reverse')
xlabel('Range (m)')
ylabel('Transmission Loss (dB)')
title(['Double Duct TL (movmean: ' num2str(k) ' datapoints) Zr=166m'])
% title(['Single Duct TL (movmean: ' num2str(k) ' datapoints) Zr=166m'])
xlim([min(min(R)) max(max(R))])
ylim([min(min(TL)) max(max(TL))])
% change legend by hand
% legend('50 Hz','300 Hz', '500 Hz','1000 Hz', '1500 Hz')   %for double
legend('50 Hz','300 Hz', '500 Hz','1000 Hz', '1500 Hz')     % for single

%%
load('combo_40_60corr.mat')
    maxcorr_lat50 = maxcorr_lat;
    maxcorr_lon50 = maxcorr_lon;
    maxcorr_val50 = cmax;
    dist50 = dist;
    size50=maxcorr_val50*200;

    load('combo_250_350corr.mat')
    maxcorr_lat300 = maxcorr_lat;
    maxcorr_lon300 = maxcorr_lon;
    maxcorr_val300 = cmax;
    dist300 = dist;
    size300=maxcorr_val300*200;

    load('combo_450_550corr.mat')
    maxcorr_lat500 = maxcorr_lat;
    maxcorr_lon500 = maxcorr_lon;
    maxcorr_val500 = cmax;
    dist500 = dist;
    size500=maxcorr_val500*200;

    load('combo_900_1100corr.mat')
    maxcorr_lat1000 = maxcorr_lat;
    maxcorr_lon1000 = maxcorr_lon;
    maxcorr_val1000 = cmax;
    dist1000 = dist;
    size1000=maxcorr_val1000*200;

    load('combo_1250_1750corr.mat')
    maxcorr_lat1500 = maxcorr_lat;
    maxcorr_lon1500 = maxcorr_lon;
    maxcorr_val1500 = cmax;
    dist1500 = dist;
    size1500=maxcorr_val1500*200;

%% SL plot
% make a plot that compares SL-TL = received
% SLs from 110-180 dB
figure
SL = [110 140 180 200];
% ind_TL=find(t >= t_beg_num & t <= t_end_num)
% dist_vec = [1500];
% for i = 1:4
% test_dist = ['dist' num2str(dist_vec(i))];
[maxdist,ind_TL] = max(dist1500);
[min_diff,ind_mindiff]=min(abs(R'-maxdist));
RL = SL-TL(ind_mindiff);
scatter(SL,RL)
hline(55)
xlabel('Source Level (dB)')
ylabel('Received Level (dB)')
title('Source Level and Received Level')
% hold on
% RL = 

