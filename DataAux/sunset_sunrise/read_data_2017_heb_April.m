close all
clear all
clc



%% Import data sunrise/sunset
filename='2017_sunrise_sunset_UTCm9_data_only';
toto=load(filename);
toto(:,1)=[];

sunrise=[];
sunset=[];
for mm=1:2:23
    sunrise=[sunrise ; toto(:,mm)];
    sunset=[sunset ; toto(:,mm+1)];
end
   



%% Data with 0 are days that don't exist (eg Feb 30)
sunrise(sunrise==9999)=[];
sunset(sunset==9999)=[];


%% Make string vectors
sunrise=num2str(sunrise, '%04u');
sunset=num2str(sunset, '%04u');

t0=datenum([2017,1,1,0,0,0]);
for tt=1:length(sunrise)
    time_sun_mat(tt)=t0+tt-1;    
end

time_sun_str=datestr(time_sun_mat, 'yyyymmdd');
time_sun_str(end,:);

%% Restric to Feb - April (data not good after April)

toto=32:32+28+31+30-1;

time_sun_mat=time_sun_mat(toto);
time_sun_str=time_sun_str(toto, :);
sunrise=sunrise(toto,:);
sunset=sunset(toto,:);

for dd=1:length(toto)
    sunrise_str(dd,:)=[time_sun_str(dd,:) 'T' sunrise(dd,:) '00'];
    sunrise_num(dd)=datenum(sunrise_str(dd,:),'yyyymmddTHHMMSS');
    sunset_str(dd,:)=[time_sun_str(dd,:) 'T' sunset(dd,:) '00'];
    sunset_num(dd)=datenum(sunset_str(dd,:),'yyyymmddTHHMMSS');
end

save sunrise_sunset_2017_feb_april sunrise_str sunrise_num sunset_str sunset_num