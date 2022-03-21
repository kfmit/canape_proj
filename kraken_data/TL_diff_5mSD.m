%% kraken trend analysis

clear all
close all
clc

addpath('/home/kfung/Downloads/CANAPE/kraken_data/5mSD_0_900km_166RD/');
fig50=open('TL_50Hz.fig');
axObjs50 = fig50.Children;
dataObjs = axObjs50.Children;
x50 = dataObjs(1).XData;
y50 = dataObjs(1).YData;

fig300=open('TL_300Hz.fig');
axObjs300 = fig300.Children;
dataObjs = axObjs300.Children;
x300 = dataObjs(1).XData;
y300 = dataObjs(1).YData;

fig500=open('TL_500Hz.fig')
axObjs500 = fig500.Children;
dataObjs = axObjs500.Children;
x500 = dataObjs(1).XData;
y500 = dataObjs(1).YData;

fig1000=open('TL_1000Hz.fig');
axObjs1000 = fig1000.Children;
dataObjs = axObjs1000.Children;
x1000 = dataObjs(1).XData;
y1000 = dataObjs(1).YData;

fig1500=open('TL_1500Hz.fig');
axObjs1500 = fig1500.Children;
dataObjs = axObjs1500.Children;
x1500 = dataObjs(1).XData;
y1500 = dataObjs(1).YData;

%%
k = 500;
R50 = movmean(x50,k);
TL50 = movmean(y50,k);

R300 = movmean(x300,k);
TL300 = movmean(y300,k);

R500 = movmean(x500,k);
TL500 = movmean(y500,k);

R1000 = movmean(x1000,k);
TL1000 = movmean(y1000,k);

R1500 = movmean(x1500,k);
TL1500 = movmean(y1500,k);

figure
hold on
set(gca,'YDir','reverse')
plot(R50,TL50)
plot(R300,TL300)
plot(R500,TL500)
plot(R1000,TL1000)
plot(R1500,TL1500)
xlabel('Range (m)')
ylabel('Transmission Loss (dB)')
title(['Moving mean, ' num2str(k) ' datapoints of TL in range, Zr=166m'])
legend('50 Hz','300 Hz', '500 Hz','1000 Hz', '1500 Hz')
