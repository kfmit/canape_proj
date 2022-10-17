%% This is a code that will plot four histograms on the same plot

close all
clear all
clc

addpath('/home/kfung/Downloads/CANAPE/mat_files/')

% change load names as needed
% load ANL_SHRU5.mat
load ANL_SHRU5_bigfreqs.mat
load sunrise_sunset_2017_feb_april.mat

% f1=[40 450 900 1250 250];
% f2=[60 550 1100 1750 350];
% 40 is EMPTY

% the four freqs i wanna look at [250 500 1000 1500]
% f1=[40 450 900 1250 250];
% f2=[60 550 1100 1750 350];
fourfreq = [5 10 15 20 25 30];
figure
tiledlayout(2,3)

for i=1:6
ff=fourfreq(i); % this picks the freq

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

avg_freq = (f1(ff)+f2(ff))/2;

nexttile

histogram(ANL_duct, ANL_vec,'Normalization', 'pdf')
hold on
histogram(ANL_no_duct, ANL_vec,'Normalization', 'pdf')
hold on
histogram(ANL_no_ice, ANL_vec,'Normalization', 'pdf')
title(['[' num2str(f1(ff)) ' - ' num2str(f2(ff)) '] Hz'])
grid on
xlim([52 100])
ylim([0 0.25])
xlabel(['ANL_{' num2str(avg_freq) '} (dB re 1 \muPa^2 / Hz)'])

% NameFig=['./new_figs/jasa_plot_hist_vs_ice_frac_new/figs_for_gif3/' num2str(f1(ff)) '-' num2str(f2(ff)) 'Hz'];
% print(gcf,NameFig,'-dpng')
    
% close(gcf)

vecsize = ANL_vec(2)-ANL_vec(1)
end

leg = legend('Ice with duct', 'Ice without duct', 'No ice')
leg.Layout.Tile = 'North'
sgtitle('Ambient Noise for Selected Frequencies')

%% another figure 500 only
% just 500

figure
for i=2:2
ff=fourfreq(i); % this picks the freq

%
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



% Histograms

% change ff as needed
% ff=5;

ANL_ice=SPL_ANL(ind_ice,ff);
ANL_no_ice=SPL_ANL(ind_no_ice,ff);

ANL_duct=ANL_ice(1:end/2);
ANL_no_duct=ANL_ice(end/2+1:end);

ANL_min=min([min(ANL_ice) min(ANL_no_ice)]);
ANL_max=max([max(ANL_ice) max(ANL_no_ice)]);
ANL_vec=linspace(ANL_min, ANL_max, 50);

avg_freq = (f1(ff)+f2(ff))/2;


histogram(ANL_duct, ANL_vec,'Normalization', 'pdf')
hold on
histogram(ANL_no_duct, ANL_vec,'Normalization', 'pdf')
hold on
histogram(ANL_no_ice, ANL_vec,'Normalization', 'pdf')
title(['Ambient Noise Level [' num2str(f1(ff)) ' - ' num2str(f2(ff)) '] Hz'])
grid on
xlim([52 92])
ylim([0 0.15])
xlabel(['ANL_{' num2str(avg_freq) '} (dB re 1 \muPa^2 / Hz)'])

% NameFig=['./new_figs/jasa_plot_hist_vs_ice_frac_new/figs_for_gif3/' num2str(f1(ff)) '-' num2str(f2(ff)) 'Hz'];
% print(gcf,NameFig,'-dpng')
    
% close(gcf)

vecsize = ANL_vec(2)-ANL_vec(1);
end

leg = legend('Ice with duct', 'Ice without duct', 'No ice')
leg.Location = 'north'