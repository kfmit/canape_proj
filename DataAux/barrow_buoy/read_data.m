close all
clear all
clc

A=load('noaa_no_header.txt');

YY=A(:,1);
MM=A(:,2);
DD=A(:,3);
hh=A(:,4);
mm=A(:,5);
temp=A(:,14);

clear A

temp(temp==999)=NaN;




for tt=1:length(YY)
    time_num(tt)=datenum(YY(tt),MM(tt),DD(tt),hh(tt),mm(tt),0);
    
end

figure
plot(time_num, temp)
datetick2('x', 'mmm.dd HH:MM')
grid on