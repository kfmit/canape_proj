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