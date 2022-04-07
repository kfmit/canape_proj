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