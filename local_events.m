clear all
close all
clc

clear all
close all
clc

load ANL_SHRU5.mat
t=timestamp_num_spectro;
SPL_ANL_nobeamform=SPL_ANL;

load ANL_SHRU5_beamform.mat

t_ice=timestamp_num_ssmi;
ice=T_ssmi.icefrac;

%% Limit to a single frequency band 250-350 Hz
ff=5;
SPL_ANL_ok=SPL_ANL_nobeamform(:,ff);
ff=2;
ANL_ver_ok=squeeze(SPL_ANL(:,ff,1));
ANL_hor_ok=squeeze(SPL_ANL(:,ff,3));

clear SPL_ANL SPL_ANL_raw SPL_ANL_nobeamform


figure
plot(t, ANL_ver_ok)
hold on
plot(t, ANL_hor_ok)
grid on
datetick('x')


%% Detect local event
s0=12; %% difference between main lobe and side lobe when beamforming vert/hor
s1=s0-3;  %% threshold in dB
s2=s0+3;
toto=find(ANL_hor_ok-ANL_ver_ok>s2 | ANL_hor_ok-ANL_ver_ok<s1);  %%% bad definition

local_ver=zeros(size(ANL_ver_ok))+NaN;
local_ver(toto)=ANL_ver_ok(toto);

local_hor=zeros(size(ANL_hor_ok))+NaN;
local_hor(toto)=ANL_hor_ok(toto);

hold on
plot(t,local_ver,'o')
hold on
plot(t,local_hor,'o')

far_hor=ANL_hor_ok;
far_hor(toto)=NaN;

far_ver=ANL_ver_ok;
far_ver(toto)=NaN;

hold on
plot(t,far_ver,'*')
hold on
plot(t,far_hor,'*')


legend('Vertical (up)', 'Horizontal','Vertical/local', 'Horizontal/local','Vertical/far', 'Horizontal/far' )

dt_anl=t(2)-t(1);
nave=round(1/dt_anl); %%% one day

figure
subplot(211)
plot(t,ANL_hor_ok-ANL_ver_ok)
grid on
datetick('x')
subplot(212)
plot(t,movmean(ANL_hor_ok-ANL_ver_ok,nave))
datetick('x')
grid on
% return

save ANL_far t far_hor far_ver
%% Number of local events per unit time
t_deb=t(1);
t_end=t(end);

% dt=0.125; %%% 3 hours (in day)
dt=1; %%% 1 day

t_int=t_deb:dt:t_end;
Nt_int=length(t_int);

n_local=zeros(size(t_int));
for tt=1:Nt_int-1
    clear toto
    toto=find(t>=t_int(tt) & t<=t_int(tt+1));
    n_local(tt)=length(find(~isnan(local_ver(toto))));
end
    
figure
plot(t_int, n_local)
grid on
datetick('x')

%%
figure

p1=subplot(611);
plot(t_int, n_local)
grid on
ylabel('Number of local event per 3 hours')
datetick('x')

p2=subplot(612);
plot(timestamp_num_ssmi,T_ssmi.icefrac)
ylabel('Ice Fraction (%)')
datetick('x')
grid on

p3=subplot(613);
plot(timestamp_num_ecmwf, T_ecmwf.W10)
ylabel('Wind speed (m/s)')
datetick('x')
grid on

p4=subplot(614);
plot(timestamp_num_ecmwf, T_ecmwf.tp)
ylabel('Total precipitation (m) ')
datetick('x')
grid on

p5=subplot(615);
plot(timestamp_num_smos, T_smos.sea_ice_thickness)
ylabel('Sea ice thickness (m) ')
datetick('x')
grid on

p6=subplot(616);
plot(timestamp_num_temp, T_temp.T2062, timestamp_num_temp, T_temp.T2061, timestamp_num_temp, T_temp.T2075)
ylabel('Underwater temperature (degree C) ')
datetick('x')
grid on

linkaxes([p1,p2,p3,p4,p5,p6],'x')
xlim([timestamp_num_spectro(1) timestamp_num_spectro(end)])