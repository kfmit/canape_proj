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

ANL_duct_avg=zeros(1,39);
ANL_no_duct_avg=zeros(1,39);
ANL_no_ice_avg=zeros(1,39);


 for i=1:39
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



%%

ANL_ice=SPL_ANL(ind_ice,ff);
ANL_no_ice=SPL_ANL(ind_no_ice,ff);

ANL_duct=ANL_ice(1:end/2);
ANL_no_duct=ANL_ice(end/2+1:end);

ANL_min=min([min(ANL_ice) min(ANL_no_ice)]);
ANL_max=max([max(ANL_ice) max(ANL_no_ice)]);
ANL_vec=linspace(ANL_min, ANL_max, 50);

avg_freq(i) = (f1(ff)+f2(ff))/2;
ANL_duct_avg(i) = mean(ANL_duct);
ANL_no_duct_avg(i) = mean(ANL_no_duct);
ANL_no_ice_avg(i)= mean(ANL_no_ice);

 end


%% Figure Creation: Trend Analysis
figure
p1 = plot(avg_freq,ANL_duct_avg,'-o','MarkerEdgeColor','#4DBEEE','MarkerFaceColor','#4DBEEE')
hold on
p2 = plot(avg_freq,ANL_no_duct_avg,'-o','MarkerEdgeColor','#D95319','MarkerFaceColor','#D95319')
p3 = plot(avg_freq,ANL_no_ice_avg,'-o','MarkerEdgeColor','#EDB120','MarkerFaceColor','#EDB120')

ylabel(['ANL (dB re 1 \muPa^2 / Hz)'])
xlabel('Frequency (Hz)')
% xticks([50 100 150 200 250 300 350 400 450 500 550 600 ...
%     650 700 750 800 850 900 950 1000 1050 1100 1150 1200 1250 1300 1350 1400 1450 ...
%     1500 1550 1600 1650 1700 1750 1800 1850 1900 1950])
xticks([100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900])
xtickangle(-45)
title('Average ANL at each frequency')
legend('Ice with duct','Ice without duct','No Ice')

%% Differences



%% Figure Creation: Histograms



% histogram fig
% figure
% histogram(ANL_duct, ANL_vec,'Normalization', 'pdf')
% hold on
% histogram(ANL_no_duct, ANL_vec,'Normalization', 'pdf')
% hold on
% histogram(ANL_no_ice, ANL_vec,'Normalization', 'pdf')
% legend('Ice with duct', 'Ice without duct', 'No ice')
% title(['Ambient noise level in [' num2str(f1(ff)) ' - ' num2str(f2(ff)) '] Hz'])
% grid on
% xlim([60 100])
% ylim([0 0.18])
% xlabel('ANL_{300} (dB re 1 \muPa^2 / Hz)')
% 
% NameFig=['./new_figs/jasa_plot_hist_vs_ice_frac_new/figs_for_gif/' num2str(f1(ff)) '-' num2str(f2(ff)) 'Hz'];
% % print(gcf,NameFig,'-dpng')
%     
% close(gcf)
% 
% % end


