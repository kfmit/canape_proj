%% Code for generating attenuation
close all
clear all
clc

c_doubleduct = [1433	1434	1438.9	1442	1444.9	1443	1441	1440.5	1442.5	1444.5	1447	1451.5	1453.9	1456.5	1458.5	1461	1462.5];
c_singleduct = [1433	1434	1440.5	1442.5	1444.5	1447	1451.5	1453.9	1456.5	1458.5	1461	1462.5];
type = 2;
f_array = [50 300 500 1000 1500];

% function [alpha1,alpha2,alpha3] = attenuation(f,c,type)
% specify unit as 1 = dB/m, 2 = dB/lambda 
% [a1,a2,a3]=attenuation(f_array(5),c_singleduct,2);

% lets use the LOWEST attn value bc we want it to work
a300 = (0.01/1000).*c_doubleduct./f_array(2);
a500 = (0.023/1000).*c_doubleduct./f_array(3);
a1000 = (0.05/1000).*c_doubleduct./f_array(4);
a1500 = (0.07/1000).*c_doubleduct./f_array(5);

%%
f_all = [50:50:2000];
% all in dB/km
alpha1=[0 0.001 0.003 0.005 0.008 0.01 0.013 0.017 0.02 0.023 0.026 0.029 0.032 0.035 0.038 0.04 0.043 ...
    0.045 0.048 0.05 0.052 0.054 0.056 0.058 0.06 0.062 0.064 0.066 0.068 0.07 0.072 0.074 0.076 0.078 ...
    0.08 0.081 0.083 0.085 0.087 0.09];  % Fisher and Simmons 77
alpha2=[0 0.002 0.003 0.006 0.009 0.012 0.016 0.02 0.024 0.028 0.032 0.036 0.04 0.043 0.047 0.051 0.054 ...
    0.057 0.061 0.064 0.067 0.07 0.073 0.076 0.078 0.081 0.084 0.087 0.09 0.092 0.095 0.098 0.101 0.103 ...
    0.106 0.109 0.112 0.115 0.118 0.121];  % Francois and Garriso 82
alpha3=[0 0.002 0.003 0.006 0.009 0.012 0.016 0.02 0.024 0.028 0.032 0.036 0.04 0.043 0.047 0.051 0.054 ...
    0.057 0.061 0.064 0.067 0.07 0.073 0.076 0.079 0.082 0.084 0.087 0.09 0.093 0.096 0.098 0.101 0.104 ...
    0.107 0.11 0.113 0.116 0.119 0.122];  % Ainslie and McColm 98
% stopped at 1.2
c_all = c_doubleduct(1);
for i=1:length(f_all)
alpha1_lam(i)= (alpha1(i)/1000)*c_all/f_all(i);
alpha2_lam(i)= (alpha2(i)/1000)*c_all/f_all(i);
alpha3_lam(i)= (alpha3(i)/1000)*c_all/f_all(i);
end

figure
plot(f_all,alpha1_lam)
hold on
plot(f_all,alpha2_lam)
plot(f_all,alpha3_lam)
xlabel('Frequency (Hz)')
ylabel('\alpha dB/\lambda')
title('Attenuation at c=1433 m/s')
legend('Fisher/Simmons','Francois/Garrison','Ainslie/McColm','Location','best')