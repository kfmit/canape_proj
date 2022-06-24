%% Code for more analysis
%% hist divider

clear all
close all
clc

addpath('/home/kfung/Downloads/CANAPE/mat_files/')
load ANLs_50_1900Hz.mat
avg_freq = avg_freq(1:38);

for i = 1:38
ANL_duct_avg(i) = mean(ANL_duct(:,i));
ANL_no_duct_avg(i) = mean(ANL_no_duct(:,i));
ANL_no_ice_avg(i)= mean(ANL_no_ice(:,i));
end

figure
p1 = plot(avg_freq(1:38),ANL_duct_avg(1:38),'-o','MarkerEdgeColor','#4DBEEE','MarkerFaceColor','#4DBEEE')
hold on
p2 = plot(avg_freq(1:38),ANL_no_duct_avg(1:38),'-o','MarkerEdgeColor','#D95319','MarkerFaceColor','#D95319')
p3 = plot(avg_freq(1:38),ANL_no_ice_avg(1:38),'-o','MarkerEdgeColor','#EDB120','MarkerFaceColor','#EDB120')

grid on
ylabel(['ANL (dB re 1 \muPa^2 / Hz)'])
xlabel('Frequency (Hz)')
% xticks([50 100 150 200 250 300 350 400 450 500 550 600 ...
%     650 700 750 800 850 900 950 1000 1050 1100 1150 1200 1250 1300 1350 1400 1450 ...
%     1500 1550 1600 1650 1700 1750 1800 1850 1900 1950])
xticks([100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900])
xtickangle(-45)
ylim([50 90])
title('Average ANL at each frequency')
legend('Ice with duct','Ice without duct','No Ice')
%% separate by MAXIMUM %%%%%%%%

ANL_duct_max = zeros(1,38);
ANL_no_duct_max = zeros(1,38); 
ANL_no_ice_max = zeros(1,38);

for i = 1:38 %used to be 38 from 50-1950
ANL_no_duct_max(i) = max(ANL_no_duct(:,i));
ANL_no_ice_max(i) = max(ANL_no_ice(:,i));
ANL_duct_max(i) = max(ANL_duct(:,i));

end

figure
p1 = plot(avg_freq,ANL_duct_max,'-o','MarkerEdgeColor','#4DBEEE','MarkerFaceColor','#4DBEEE')
hold on
p2 = plot(avg_freq,ANL_no_duct_max,'-o','MarkerEdgeColor','#D95319','MarkerFaceColor','#D95319')
p3 = plot(avg_freq,ANL_no_ice_max,'-o','MarkerEdgeColor','#EDB120','MarkerFaceColor','#EDB120')

ylabel(['Max ANL (dB re 1 \muPa^2 / Hz)'])
xlabel('Frequency (Hz)')
grid on
% xticks([50 100 150 200 250 300 350 400 450 500 550 600 ...
%     650 700 750 800 850 900 950 1000 1050 1100 1150 1200 1250 1300 1350 1400 1450 ...
%     1500 1550 1600 1650 1700 1750 1800 1850 1900 1950])
xticks([100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900])
xtickangle(-45)
title('Max ANL at each frequency')
legend('Ice with duct','Ice without duct','No Ice')
%% separate by min

ANL_duct_min = zeros(1,38);
ANL_no_duct_min = zeros(1,38); 
ANL_no_ice_min = zeros(1,38);
for i = 1:38
ANL_no_duct_min(i) = min(ANL_no_duct(:,i));
ANL_no_ice_min(i) = min(ANL_no_ice(:,i));
ANL_duct_min(i) = min(ANL_duct(:,i));

end

figure
p1 = plot(avg_freq,ANL_duct_min,'-o','MarkerEdgeColor','#4DBEEE','MarkerFaceColor','#4DBEEE')
hold on
p2 = plot(avg_freq,ANL_no_duct_min,'-o','MarkerEdgeColor','#D95319','MarkerFaceColor','#D95319')
p3 = plot(avg_freq,ANL_no_ice_min,'-o','MarkerEdgeColor','#EDB120','MarkerFaceColor','#EDB120')

ylabel(['Min ANL (dB re 1 \muPa^2 / Hz)'])
xlabel('Frequency (Hz)')
grid on
% xticks([50 100 150 200 250 300 350 400 450 500 550 600 ...
%     650 700 750 800 850 900 950 1000 1050 1100 1150 1200 1250 1300 1350 1400 1450 ...
%     1500 1550 1600 1650 1700 1750 1800 1850 1900 1950])
xticks([100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900])
xtickangle(-45)
title('Min ANL at each frequency')
legend('Ice with duct','Ice without duct','No Ice')

%% Max and min plot bc why not

figure
p1 = plot(avg_freq,ANL_duct_max-ANL_duct_min,'-o','MarkerEdgeColor','#4DBEEE','MarkerFaceColor','#4DBEEE')
hold on
p2 = plot(avg_freq,ANL_no_duct_max-ANL_no_duct_min,'-o','MarkerEdgeColor','#D95319','MarkerFaceColor','#D95319')
p3 = plot(avg_freq,ANL_no_ice_max-ANL_no_ice_min,'-o','MarkerEdgeColor','#EDB120','MarkerFaceColor','#EDB120')

ylabel(['Min ANL (dB re 1 \muPa^2 / Hz)'])
xlabel('Frequency (Hz)')
grid on
% xticks([50 100 150 200 250 300 350 400 450 500 550 600 ...
%     650 700 750 800 850 900 950 1000 1050 1100 1150 1200 1250 1300 1350 1400 1450 ...
%     1500 1550 1600 1650 1700 1750 1800 1850 1900 1950])
xticks([100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900])
xtickangle(-45)
title('Difference of Max-Min ANL at each frequency')
legend('Ice with duct','Ice without duct','No Ice')

%% separate by median
ANL_duct_med = zeros(1,38);
ANL_no_duct_med = zeros(1,38); 
ANL_no_ice_med = zeros(1,38);
for i = 1:38
ANL_no_duct_med(i) = median(ANL_no_duct(:,i));
ANL_no_ice_med(i) = median(ANL_no_ice(:,i));
ANL_duct_med(i) = median(ANL_duct(:,i));
end


figure
p1 = plot(avg_freq,ANL_duct_med,'-o','MarkerEdgeColor','#4DBEEE','MarkerFaceColor','#4DBEEE')
hold on
p2 = plot(avg_freq,ANL_no_duct_med,'-o','MarkerEdgeColor','#D95319','MarkerFaceColor','#D95319')
p3 = plot(avg_freq,ANL_no_ice_med,'-o','MarkerEdgeColor','#EDB120','MarkerFaceColor','#EDB120')
grid on
ylabel(['Median of ANL (dB re 1 \muPa^2 / Hz)'])
xlabel('Frequency (Hz)')
% xticks([50 100 150 200 250 300 350 400 450 500 550 600 ...
%     650 700 750 800 850 900 950 1000 1050 1100 1150 1200 1250 1300 1350 1400 1450 ...
%     1500 1550 1600 1650 1700 1750 1800 1850 1900 1950])
xticks([100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900])
xtickangle(-45)
title('Median of ANL at each frequency')
legend('Ice with duct','Ice without duct','No Ice')




%% separate by biggest bin (this is mode right) - yes
% Total Variation Distance can be put into this too

% INITIALIZE
max_duct=zeros(1,38);
max_no_duct=zeros(1,38);
max_no_ice=zeros(1,38);
ind_duct = zeros(1,38);
ind_no_duct = zeros(1,38);
ind_no_ice = zeros(1,38);
totalvar_dist_duct_noduct = zeros(1,38);
totalvar_dist_duct_noice = zeros(1,38);
totalvar_dist_noice_noduct = zeros(1,38);

ANL_min=min([min(ANL_ice) min(ANL_no_ice)]);
ANL_max=max([max(ANL_ice) max(ANL_no_ice)]);
ANL_vec=linspace(ANL_min, ANL_max, 50);

% histrograms with 50 bins, N is the count in each bin and thers also edges
% speicfy bin vector

for i=1:38
    % number in each bin, edges
[N_duct,edges_duct] = histcounts(ANL_duct(:,i),ANL_vec,'Normalization','pdf');    %,'Normalization','pdf');
[N_no_duct,edges_no_duct] = histcounts(ANL_no_duct(:,i),ANL_vec,'Normalization','pdf');    %,'Normalization','pdf');
[N_no_ice,edges_no_ice] = histcounts(ANL_no_ice(:,i),ANL_vec,'Normalization','pdf');    %,'Normalization','pdf'); %,'Normalization','pdf');

% edge realignment
edges_duct=edges_duct(1:end-1)+(1.369/2);
edges_no_duct=edges_no_duct(1:end-1)+(1.369/2);
edges_no_ice=edges_no_ice(1:end-1)+(1.369/2);

% index refers to which bin of ANL_vec its in
% max_ the largest bing, ind_duct = index of largest bin
[max_duct(i), ind_duct(i)]=max(N_duct);
[max_no_duct(i), ind_no_duct(i)]=max(N_no_duct);
[max_no_ice(i), ind_no_ice(i)]=max(N_no_ice);

% total var difference
% N_* is rewritten everytime and represents the count
totalvar_dist_duct_noduct(i)=sum(abs(N_duct-N_no_duct))/2;
totalvar_dist_duct_noice(i) = sum(abs(N_duct-N_no_ice))/2;
totalvar_dist_noice_noduct(i) = sum(abs(N_no_ice-N_no_duct))/2;

% pairwise diff of max modes and screw preallocating
pair_dist_duct_noduct(i)=edges_duct(ind_duct(i))-edges_no_duct(ind_no_duct(i));
pair_dist_duct_noice(i)=edges_no_ice(ind_no_ice(i))-edges_duct(ind_duct(i));
pair_dist_noice_noduct(i)=edges_no_ice(ind_no_ice(i))-edges_no_duct(ind_no_duct(i));
end

%% Figure of Total Variation Dist %%%%%%%%%%%%%%%%%
figure
p1 = plot(avg_freq,totalvar_dist_duct_noduct,'-o','Color','#7E2F8E','MarkerEdgeColor','#7E2F8E','MarkerFaceColor','#7E2F8E');
hold on
p2 = plot(avg_freq,totalvar_dist_duct_noice,'-o','Color','#77AC30','MarkerEdgeColor','#77AC30','MarkerFaceColor','#77AC30');
p3 = plot(avg_freq,totalvar_dist_noice_noduct,'-o','Color','#FF8800','MarkerEdgeColor','#FF8800','MarkerFaceColor','#FF8800');

grid on
ylabel(['Total variation Distance'])
xlabel('Frequency (Hz)')
xticks([100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900])
xtickangle(-45)
title('Total Variation Distance of ANL at each frequency')
legend('Ice with duct/Ice without Duct','Ice with duct/No Ice','Ice without Duct/No Ice','Location','best')

%% MODE DIFFS - this is essentially pairwise differnce %%%%%%%%%%%%%%%%%%%%%%
figure
p1 = plot(avg_freq,ANL_vec(ind_duct)-ANL_vec(ind_no_duct),'-o','Color','#7E2F8E','MarkerEdgeColor','#7E2F8E','MarkerFaceColor','#7E2F8E');
hold on
p2 = plot(avg_freq,ANL_vec(ind_no_ice)-ANL_vec(ind_duct),'-o','Color','#77AC30','MarkerEdgeColor','#77AC30','MarkerFaceColor','#77AC30');
p3 = plot(avg_freq,ANL_vec(ind_no_ice)-ANL_vec(ind_no_duct),'-o','Color','#FF8800','MarkerEdgeColor','#FF8800','MarkerFaceColor','#FF8800');

grid on
ylabel(['Pairwise Differnce of ANL (dB re 1 \muPa^2 / Hz)'])
xlabel('Frequency (Hz)')
xticks([100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900])
xtickangle(-45)
ylim([-2 25])
title('Pairwise Difference of ANL between Modes of each frequency')
legend('Ice with duct/Ice without duct','Ice with duct/No Ice','Ice without duct/No Ice','Location','best')

%% MODES modes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
p1 = plot(avg_freq,ANL_vec(ind_duct),'-o','MarkerEdgeColor','#4DBEEE','MarkerFaceColor','#4DBEEE');
hold on
p2 = plot(avg_freq,ANL_vec(ind_no_duct),'-o','MarkerEdgeColor','#D95319','MarkerFaceColor','#D95319');
p3 = plot(avg_freq,ANL_vec(ind_no_ice),'-o','MarkerEdgeColor','#EDB120','MarkerFaceColor','#EDB120');

grid on
ylabel(['ANL dB re 1 \muPa^2 / Hz'])
xlabel('Frequency (Hz)')
xticks([100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900])
xtickangle(-45)
ylim([50 90])
title('Peak Probability of ANL at each frequency')
legend('Ice with duct','Ice without duct','No Ice')

%% LOG MODE PLOTS %%%%%%%%%%%
figure
p1 = plot(avg_freq,log10(ANL_vec(ind_duct)),'-o','MarkerEdgeColor','#4DBEEE','MarkerFaceColor','#4DBEEE');
hold on
p2 = plot(avg_freq,log10(ANL_vec(ind_no_duct)),'-o','MarkerEdgeColor','#D95319','MarkerFaceColor','#D95319');
p3 = plot(avg_freq,log10(ANL_vec(ind_no_ice)),'-o','MarkerEdgeColor','#EDB120','MarkerFaceColor','#EDB120');

grid on
ylabel(['Mode of ANL log(dB re 1 \muPa^2 / Hz)'])
xlabel('Frequency (Hz)')
xticks([100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900])
xtickangle(-45)
title('Log Mode of ANL at each frequency')
legend('Ice with duct','Ice without duct','No Ice')

%% Try by pairwise again PAIRWISE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
p1 = plot(avg_freq,pair_dist_duct_noduct,'-o','Color','#7E2F8E','MarkerEdgeColor','#7E2F8E','MarkerFaceColor','#7E2F8E');
hold on
p2 = plot(avg_freq,pair_dist_duct_noice,'-o','Color','#77AC30','MarkerEdgeColor','#77AC30','MarkerFaceColor','#77AC30');
p3 = plot(avg_freq,pair_dist_noice_noduct,'-o','Color','#FF8800','MarkerEdgeColor','#FF8800','MarkerFaceColor','#FF8800');

grid on
ylabel(['Pairwise Differnce of Peak (dB re 1 \muPa^2 / Hz)'])
xlabel('Frequency (Hz)')
ylim([-2 25])
xticks([100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900])
xtickangle(-45)
title('Pairwise Difference of ANL at each frequency')
legend('Ice with duct/Ice without Duct','Ice with Duct/No Ice','Ice without Duct/No Ice','Location','best')

%% Pairwise but just the last on %%%%%%%%%%%%%%%%%%%%%%%e
figure
p1 = plot(avg_freq,pair_dist_duct_noduct,'-o','Color','#FF8800','MarkerEdgeColor','#FF0000','MarkerFaceColor','#FF0000');

ylabel(['Pairwise Differnce of Mode (dB re 1 \muPa^2 / Hz)'])
xlabel('Frequency (Hz)')
xticks([100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900])
xtickangle(-45)
title('Pairwise Difference of ANL at each frequency')
legend('Ice with duct/Ice without Duct','Location','best')