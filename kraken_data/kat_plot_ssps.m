%% ssps for me
clc
clear all
close all

% should i include or not include 3000
double_depth = [0 40 50 55 70 85 127 160 210 226 250 277 317 353 420 520 780];
double_ss = [1433 1434 1438.9 1442 1444.9 1443 1441 1440.5 1442.5 1444.5 1447 1451.5 1453.9 1456.5 1458.5 1461 1462.5];
single_depth = [0 40 160 210 226 250 277 317 353 420 520 780];
single_ss = [1433 1434 1440.5 1442.5 1444.5 1447 1451.5 1453.9 1456.5 1458.5 1461 1462.5];

figure
tiledlayout(1,2)
nexttile
plot(double_ss,double_depth)
xlabel('Sound Speed (m/s)')
ylabel('Depth (m)')
set(gca,'YDir', 'reverse')
title('Double Duct SSP')

nexttile
plot(single_ss,single_depth)
xlabel('Sound Speed (m/s)')
ylabel('Depth (m)')
set(gca,'YDir', 'reverse')
title('Single Duct SSP')