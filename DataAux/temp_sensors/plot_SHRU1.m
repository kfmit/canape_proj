close all
clear all
clc

name_shru='SHRU1';
tpod={'T2085', 'T2073', 'T2050'};
sbe={'sbe3125'};

temp_shru=NaN(5e6,4);
tm_shru=NaN(5e6,4);

for tt=1:length(tpod)

    tpod_file=['./tpods/' tpod{tt} '.out'];

    tdata = load(tpod_file);

    day = tdata(:,2);
    mon =  tdata(:,3);
    yr =  2000+tdata(:,4);
    hr = tdata(:,5);
    mins = tdata(:,6);
    sec = tdata(:,7);

    temp = tdata(:,8);
    tm = datenum(yr,mon,day,hr,mins,sec);
    

    temp_shru(1:length(temp),tt)=temp;
    tm_shru(1:length(tm),tt)=tm;
%     plot(tm,temp) 
%     grid on
%     datetick('x','keeplimits');
%     xlabel('Date');
%     ylabel('Temperature (deg C)')

end
tt=tt+1;


sbe_file=['./SBE/' sbe{1} '.dat'];

D = load(sbe_file);

temp = D(:,1);
P = D(:,2);
day = D(:,3);
mon = D(:,4);
yr = D(:,5);
hr = D(:,6);
mm = D(:,7);
ss= D(:,8);

tm = datenum(yr,mon,day,hr,mm,ss);


% plot(tm,temp) 
% grid on
% datetick('x','keeplimits');
% xlabel('Date');
% ylabel('Temperature (deg C)')

temp_shru(1:length(temp),tt)=temp;
tm_shru(1:length(tm),tt)=tm;


figure
plot(tm_shru,temp_shru)
datetick('x','keeplimits');
xlabel('Date');
ylabel('Temperature (deg C)')
grid on

%% Put everybody on the same basis

t_start_vec=[2016,10,22,0,0,0];
t_end_vec  =[2017,09,30,0,0,0];

tm_start=datenum(t_start_vec);
tm_end=datenum(t_end_vec);

tm0=tm_start:2/60/24:tm_end;
temp_shru_ok=zeros(length(tm0),tt);

T=table(datestr(tm0,'yyyymmddHHMMSS'),'VariableNames',{'timestamp'});

for ii=1:tt-1
    toto=find(tm_shru(:,ii) > tm_start & tm_shru(:,ii) < tm_end);
    temp_shru_ok(:,ii)=interp1(tm_shru(toto,ii),temp_shru(toto,ii),tm0,'linear');  
    T=[T table(temp_shru_ok(:,ii), 'VariableNames',{tpod{ii}})];
end

ii=tt;
toto=find(tm_shru(:,ii) > tm_start & tm_shru(:,ii) < tm_end);
temp_shru_ok(:,ii)=interp1(tm_shru(toto,ii),temp_shru(toto,ii),tm0,'linear');  
T=[T table(temp_shru_ok(:,ii), 'VariableNames',sbe)];


hold on
plot(tm0,temp_shru_ok,'o')

writetable(T,['variables_temp_' name_shru '.csv']) 

%%
figure
plot(tm0,temp_shru_ok)
datetick('x','keeplimits');
xlabel('Date');
ylabel('Temperature (deg C)')
grid on