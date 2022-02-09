%% Code for more analysis
%% hist divider

clear all
close all
clc

load ANLs_50_1900Hz.mat

%% separate by max

ANL_duct_max = zeros(1,39);
ANL_no_duct_max = zeros(1,39); 
ANL_no_ice_max = zeros(1,39);

for i = 1:39
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
% xticks([50 100 150 200 250 300 350 400 450 500 550 600 ...
%     650 700 750 800 850 900 950 1000 1050 1100 1150 1200 1250 1300 1350 1400 1450 ...
%     1500 1550 1600 1650 1700 1750 1800 1850 1900 1950])
xticks([100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900])
xtickangle(-45)
title('Max ANL at each frequency')
legend('Ice with duct','Ice without duct','No Ice')
%% separate by min

ANL_duct_min = zeros(1,39);
ANL_no_duct_min = zeros(1,39); 
ANL_no_ice_min = zeros(1,39);
for i = 1:39
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
% xticks([50 100 150 200 250 300 350 400 450 500 550 600 ...
%     650 700 750 800 850 900 950 1000 1050 1100 1150 1200 1250 1300 1350 1400 1450 ...
%     1500 1550 1600 1650 1700 1750 1800 1850 1900 1950])
xticks([100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900])
xtickangle(-45)
title('Difference of Max-Min ANL at each frequency')
legend('Ice with duct','Ice without duct','No Ice')

%% separate by median
ANL_duct_med = zeros(1,39);
ANL_no_duct_med = zeros(1,39); 
ANL_no_ice_med = zeros(1,39);
for i = 1:39
ANL_no_duct_med(i) = median(ANL_no_duct(:,i));
ANL_no_ice_med(i) = median(ANL_no_ice(:,i));
ANL_duct_med(i) = median(ANL_duct(:,i));
end


figure
p1 = plot(avg_freq,ANL_duct_med,'-o','MarkerEdgeColor','#4DBEEE','MarkerFaceColor','#4DBEEE')
hold on
p2 = plot(avg_freq,ANL_no_duct_med,'-o','MarkerEdgeColor','#D95319','MarkerFaceColor','#D95319')
p3 = plot(avg_freq,ANL_no_ice_med,'-o','MarkerEdgeColor','#EDB120','MarkerFaceColor','#EDB120')

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

max_duct=zeros(1,39);
max_no_duct=zeros(1,39);
max_no_ice=zeros(1,39);
ind_duct = zeros(1,39);
ind_no_duct = zeros(1,39);
ind_no_ice = zeros(1,39);

ANL_min=min([min(ANL_ice) min(ANL_no_ice)]);
ANL_max=max([max(ANL_ice) max(ANL_no_ice)]);
ANL_vec=linspace(ANL_min, ANL_max, 50);

% histrograms with 15 bins, N is the count in each bin and thers also edges
% speicfy bin vector

for i=1:39
[N_duct,edges_duct] = histcounts(ANL_duct(:,i),ANL_vec,'Normalization','pdf');
[N_no_duct,edges_no_duct] = histcounts(ANL_no_duct(:,i),ANL_vec,'Normalization','pdf');
[N_no_ice,edges_no_ice] = histcounts(ANL_no_ice(:,i),ANL_vec,'Normalization','pdf');

% index refers to which bin of ANL_vec its in
[max_duct(i), ind_duct(i)]=max(N_duct);
[max_no_duct(i), ind_no_duct(i)]=max(N_no_duct);
[max_no_ice(i), ind_no_ice(i)]=max(N_no_ice);

end

%%
figure
p1 = plot(avg_freq,ANL_vec(ind_duct),'-o','MarkerEdgeColor','#4DBEEE','MarkerFaceColor','#4DBEEE')
hold on
p2 = plot(avg_freq,ANL_vec(ind_no_duct),'-o','MarkerEdgeColor','#D95319','MarkerFaceColor','#D95319')
p3 = plot(avg_freq,ANL_vec(ind_no_ice),'-o','MarkerEdgeColor','#EDB120','MarkerFaceColor','#EDB120')

ylabel(['Mode of ANL (dB re 1 \muPa^2 / Hz)'])
xlabel('Frequency (Hz)')
% xticks([50 100 150 200 250 300 350 400 450 500 550 600 ...
%     650 700 750 800 850 900 950 1000 1050 1100 1150 1200 1250 1300 1350 1400 1450 ...
%     1500 1550 1600 1650 1700 1750 1800 1850 1900 1950])
xticks([100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900])
xtickangle(-45)
title('Mode of ANL at each frequency')
legend('Ice with duct','Ice without duct','No Ice')




%% separate by total variation distance of probaility measures
% use the histcounts outputs to calculate

