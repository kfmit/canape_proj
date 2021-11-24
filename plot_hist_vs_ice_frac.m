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
title('Ambient noise level')
p2=subplot(212);
plot(timestamp_num_ssmi,T_ssmi.icefrac)
title('Ice concentration')
grid on
datetick('x')
linkaxes([p1 p2], 'x')

%% Histogram

for ff=1:Nf
    ANL_ice=SPL_ANL(ind_ice,ff);
    ANL_no_ice=SPL_ANL(ind_no_ice,ff);
    SPL_ice=SPL_raw(ind_ice,ff);
    SPL_no_ice=SPL_raw(ind_no_ice,ff);
    
    ANL_min=min([min(ANL_ice) min(ANL_no_ice)]);
    ANL_max=max([max(ANL_ice) max(ANL_no_ice)]);
    ANL_vec=linspace(ANL_min, ANL_max, 50);
    
    SPL_min=min([min(SPL_ice) min(SPL_no_ice)]);
    SPL_max=max([max(SPL_ice) max(SPL_no_ice)]);
    SPL_vec=linspace(SPL_min, SPL_max, 50);
    
    figure
    subplot(211)
    histogram(ANL_ice, ANL_vec,'Normalization', 'pdf')
    hold on
    histogram(ANL_no_ice, ANL_vec,'Normalization', 'pdf')
    title(['Ambient noise level in [' num2str(f1(ff)) ' - ' num2str(f2(ff)) '] Hz'])
    legend('Ice', 'No Ice')
    subplot(212)
    histogram(SPL_ice, SPL_vec,'Normalization', 'pdf')
    hold on
    histogram(SPL_no_ice, SPL_vec,'Normalization', 'pdf')
    legend('Ice', 'No Ice')
    title(['Recorded sound level in [' num2str(f1(ff)) ' - ' num2str(f2(ff)) '] Hz'])
    
end

