close all
clear all
clc

%%%% NB: local time is UTC-8


%% Import data sunrise/sunset
filename='2016_sunrise_sunset_data_only';
toto=load(filename);
toto(:,1)=[];

sunrise=[];
sunset=[];
for mm=1:2:23
    sunrise=[sunrise ; toto(:,mm)];
    sunset=[sunset ; toto(:,mm+1)];
end
   
filename='2017_sunrise_sunset_data_only';
toto=load(filename);
toto(:,1)=[];

for mm=1:2:23
    sunrise=[sunrise ; toto(:,mm)];
    sunset=[sunset ; toto(:,mm+1)];
end
   
%% Import data twilight

filename='2016_astro_twilight_data_only';
toto=load(filename);
toto(:,1)=[];

twi_beg=[];
twi_end=[];
for mm=1:2:23
    twi_beg=[twi_beg ; toto(:,mm)];
    twi_end=[twi_end ; toto(:,mm+1)];
end
 
filename='2017_astro_twilight_data_only';
toto=load(filename);
toto(:,1)=[];

for mm=1:2:23
    twi_beg=[twi_beg ; toto(:,mm)];
    twi_end=[twi_end ; toto(:,mm+1)];
end
    


%% Data with 0 are days that don't exist (eg Feb 30)
sunrise(sunrise==0)=[];
sunset(sunset==0)=[];

twi_beg(twi_beg==0)=[];
twi_end(twi_end==0)=[];

%% Make string vectors
sunrise=num2str(sunrise, '%04u');
sunset=num2str(sunset, '%04u');
twi_beg=num2str(twi_beg, '%04u');
twi_end=num2str(twi_end, '%04u');


t0=datenum([2016,1,1,0,0,0]);
for tt=1:length(sunrise)
    time_sun_mat(tt)=t0+tt-1;    
end
time_sun_str=datestr(time_sun_mat, 'yyyymmdd');

save variables_sun sunrise sunset twi_beg twi_end time_sun_str time_sun_mat
 