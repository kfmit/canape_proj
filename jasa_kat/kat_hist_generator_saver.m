%% This code creates

close all
clear all
clc

% change load names as needed
% load ANL_SHRU5.mat
load ANL_SHRU5_bigfreqs.mat
load sunrise_sunset_2017_feb_april.mat

% f1=[40 450 900 1250 250];
% f2=[60 550 1100 1750 350];
% 40 is EMPTY

for i=4:22
ff=i; % this picks the freq

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
binsize(i) = ANL_vec(2)-ANL_vec(1)

avg_freq = (f1(ff)+f2(ff))/2;

figure
histogram(ANL_duct, ANL_vec,'Normalization', 'pdf')
hold on
histogram(ANL_no_duct, ANL_vec,'Normalization', 'pdf')
hold on
histogram(ANL_no_ice, ANL_vec,'Normalization', 'pdf')
legend('Ice with duct', 'Ice without duct', 'No ice')
title(['Ambient noise level in [' num2str(f1(ff)) ' - ' num2str(f2(ff)) '] Hz'])
grid on
xlim([50 105])
ylim([0 0.20])
xlabel(['ANL_{' num2str(avg_freq) '} (dB re 1 \muPa^2 / Hz)'])

NameFig=['./new_figs/jasa_plot_hist_vs_ice_frac_new/figs_for_gif3/' num2str(f1(ff)) '-' num2str(f2(ff)) 'Hz'];
print(gcf,NameFig,'-dpng')
    
close(gcf)

end


