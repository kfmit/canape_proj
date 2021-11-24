close all
clear all
clc

load ANL_SHRU5.mat
load sunrise_sunset_2017_feb_april.mat


ff=2;

figure
subplot(211)
yyaxis left
plot(timestamp_num_spectro,SPL_ANL(:,ff))
yyaxis right
plot(timestamp_num_ssmi,T_ssmi.icefrac)
ylim([-10 110])
grid on
datetick('x')
title('Ambient noise level')

subplot(212)
yyaxis left
plot(timestamp_num_spectro,SPL_raw(:,ff))
yyaxis right
plot(timestamp_num_ssmi,T_ssmi.icefrac)
ylim([-10 110])
grid on
datetick('x','mmmyy')
title('Raw data')

%%

for ff=1:Nf

    SPL_ANL_day=[];
    SPL_ANL_night=[];
    SPL_raw_day=[];
    SPL_raw_night=[];

    %%%% acoustics is in UTC time and sunrise in UTC-9
    local_to_UTC=-9/24;


    all_day=[];
    all_night=[];
    for dd=1:length(sunset_num)-1

        sunrise_today=sunrise_num(dd)+local_to_UTC;
        sunset_today=sunset_num(dd)+local_to_UTC;
        sunrise_tomorrow=sunrise_num(dd+1)+local_to_UTC;
        dt=0/24; %%% margin before/after sunset/sunrise

        day=find(timestamp_num_spectro > sunrise_today+dt & timestamp_num_spectro < sunset_today-dt);
        night=find(timestamp_num_spectro > sunset_today +dt & timestamp_num_spectro < sunrise_tomorrow -dt);

        all_day=[all_day ; day];
        all_night=[all_night; night];

        SPL_ANL_day=[SPL_ANL_day ; SPL_ANL(day,ff)];
        SPL_ANL_night=[SPL_ANL_night ; SPL_ANL(night,ff)];

        SPL_raw_day=[SPL_raw_day ; SPL_raw(day,ff)];
        SPL_raw_night=[SPL_raw_night ; SPL_raw(night,ff)];

    end

    ANL_min=min([min(SPL_ANL_day) min(SPL_ANL_night)]);
    ANL_max=max([max(SPL_ANL_day) max(SPL_ANL_night)]);
    ANL_vec=linspace(ANL_min, ANL_max, 50);
    
    SPL_min=min([min(SPL_raw_day) min(SPL_raw_night)]);
    SPL_max=max([max(SPL_raw_day) max(SPL_raw_night)]);
    SPL_vec=linspace(SPL_min, SPL_max, 50);
    
    figure
    subplot(211)
    histogram(SPL_ANL_day, ANL_vec,'Normalization', 'pdf')
    hold on
    histogram(SPL_ANL_night, ANL_vec,'Normalization', 'pdf')
    title(['Ambient noise level in [' num2str(f1(ff)) ' - ' num2str(f2(ff)) '] Hz'])
    legend('day', 'night')
    subplot(212)
    histogram(SPL_raw_day, SPL_vec,'Normalization', 'pdf')
    hold on
    histogram(SPL_raw_night, SPL_vec,'Normalization', 'pdf')
    legend('day', 'night')
    title(['Recorded sound level in [' num2str(f1(ff)) ' - ' num2str(f2(ff)) '] Hz'])
    
end

%% Plot time series

ff=2;
figure
plot(timestamp_num_spectro, SPL_ANL(:,ff))
hold on
plot(timestamp_num_spectro(all_day), SPL_ANL(all_day,ff), '.')
grid on
xlim([sunrise_num(1) sunset_num(end)])
% datetick2('x', 'mmm-dd HH:MM')


%% Link with Tsurf?
% ff=2;
% 
% figure
% plot(timestamp_num_smos, T_smos.Tsurf)
% grid on
% % xlim([sunrise_num(1) sunrise_num(end)]+local_to_UTC)
% datetick2('x', 'mmm-dd HH:MM')
%  
% %%% restric to period of interest
% % toto=find(timestamp_num_smos > sunrise_num(1)+local_to_UTC & timestamp_num_smos < sunrise_num(end)+local_to_UTC);
% toto=1:length(timestamp_num_smos);
% Tsurf_ok=T_smos.Tsurf(toto);
% time_T=timestamp_num_smos(toto);
% 
% temp_inc=find(diff(Tsurf_ok)>0);
% temp_dec=find(diff(Tsurf_ok)<0);
% temp_dec_fast=find(diff(Tsurf_ok)<-0);
% 
% figure
% plot(time_T, Tsurf_ok)
% hold on
% plot(time_T(temp_inc), Tsurf_ok(temp_inc), '*')
% plot(time_T(temp_dec), Tsurf_ok(temp_dec), '*')
% plot(time_T(temp_dec_fast), Tsurf_ok(temp_dec_fast), 'o')
% datetick2('x', 'mmm-dd HH:MM')
% 
% t_anl=timestamp_num_spectro;
% 
% toto=[];
% for tt=1:length(temp_dec_fast)
%     toto=[toto ; find(t_anl > time_T(temp_dec_fast(tt)) & t_anl < time_T(temp_dec_fast(tt)+1))];
% end
% 
% 
% 
% figure
% p1=subplot(211);
% plot(t_anl, SPL_ANL(:,ff))
% hold on
% plot(t_anl(toto), SPL_ANL(toto,ff), '.')
% grid on
% datetick2('x', 'mmm-dd HH:MM')
% p2=subplot(212);
% plot(time_T, Tsurf_ok)
% hold on
% plot(time_T(temp_dec_fast), Tsurf_ok(temp_dec_fast), 'o')
% grid on
% datetick2('x', 'mmm-dd HH:MM')
% linkaxes([p1 p2], 'x')
% 
