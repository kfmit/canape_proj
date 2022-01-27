close all
clear all
clc

% change load names as needed
% load ANL_SHRU5.mat
load ANL_SHRU5_newfreq.mat
load sunrise_sunset_2017_feb_april.mat

% f1=[40 450 900 1250 250];
% f2=[60 550 1100 1750 350];

ff=5; % this picks the freq

figure
subplot(211)
yyaxis left
plot(timestamp_num_spectro,SPL_ANL(:,ff))
yyaxis right
plot(timestamp_num_ssmi,T_ssmi.icefrac)
ylim([-10 110])
grid on
datetick('x')
title(['Ambient noise level for [' num2str(f1(ff)) ' - ' num2str(f2(ff)) ']  Hz'])

subplot(212)
yyaxis left
plot(timestamp_num_spectro,SPL_raw(:,ff))
yyaxis right
plot(timestamp_num_ssmi,T_ssmi.icefrac)
ylim([-10 110])
grid on
datetick('x','mmmyy')
title(['Raw data for [' num2str(f1(ff)) ' - ' num2str(f2(ff)) ']  Hz'])

%%
Nice=length(timestamp_num_ssmi);
%%%%% Period with ice
toto=find(T_ssmi.icefrac>85);
beg_ice=timestamp_num_ssmi(toto(1));
end_ice=timestamp_num_ssmi(toto(end));

ind_ice=find(timestamp_num_spectro>beg_ice & timestamp_num_spectro<end_ice);

%%%%% Period without ice
toto=find(T_ssmi.icefrac(1:floor(Nice/2))<15);
end_no_ice=timestamp_num_ssmi(toto(end));
ind_no_ice_1=find(timestamp_num_spectro < end_no_ice);

toto=find(T_ssmi.icefrac(end:-1:floor(Nice/2))<15);
beg_no_ice=timestamp_num_ssmi(Nice-toto(end));
ind_no_ice_2=find(timestamp_num_spectro > beg_no_ice);

ind_no_ice=[ind_no_ice_1 ; ind_no_ice_2];


figure
p1=subplot(211);
plot(timestamp_num_spectro,SPL_ANL(:,ff))
hold on
plot(timestamp_num_spectro(ind_ice),SPL_ANL(ind_ice,ff))
plot(timestamp_num_spectro(ind_no_ice),SPL_ANL(ind_no_ice,ff), '.')
grid on
datetick('x')
title(['Ambient noise level for [' num2str(f1(ff)) ' - ' num2str(f2(ff)) ']  Hz'])
p2=subplot(212);
plot(timestamp_num_ssmi,T_ssmi.icefrac)
title('Ice concentration')
grid on
datetick('x')
linkaxes([p1 p2], 'x')



%% Histograms

% change ff as needed
% ff=5;

ANL_ice=SPL_ANL(ind_ice,ff);
ANL_no_ice=SPL_ANL(ind_no_ice,ff);

ANL_duct=ANL_ice(1:end/2);
ANL_no_duct=ANL_ice(end/2+1:end);

ANL_min=min([min(ANL_ice) min(ANL_no_ice)]);
ANL_max=max([max(ANL_ice) max(ANL_no_ice)]);
ANL_vec=linspace(ANL_min, ANL_max, 50);



figure
histogram(ANL_duct, ANL_vec,'Normalization', 'pdf')
hold on
histogram(ANL_no_duct, ANL_vec,'Normalization', 'pdf')
hold on
histogram(ANL_no_ice, ANL_vec,'Normalization', 'pdf')
legend('Ice with duct', 'Ice without duct', 'No ice')
title(['Ambient noise level in [' num2str(f1(ff)) ' - ' num2str(f2(ff)) '] Hz'])
grid on
xlim([60 105])
xlabel('ANL_{300} (dB re 1 \muPa^2 / Hz)')

