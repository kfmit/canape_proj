%% A code to make correlation and covariance plots of ANL

close all
clear all
clc

addpath('/home/kfung/Downloads/CANAPE/mat_files');

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
ANL_duct=zeros(21488,39);
ANL_no_duct=zeros(21488,39);
ANL_no_ice=zeros(7180,39);


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
ANL_no_ice(:,i)=SPL_ANL(ind_no_ice,ff);

ANL_duct(:,i)=ANL_ice(1:end/2);
ANL_no_duct(:,i)=ANL_ice(end/2+1:end);

ANL_min=min([min(ANL_ice) min(ANL_no_ice)]);
ANL_max=max([max(ANL_ice) max(ANL_no_ice)]);
ANL_vec=linspace(ANL_min, ANL_max, 50);

avg_freq(i) = (f1(ff)+f2(ff))/2;
ANL_duct_avg(i) = mean(ANL_duct(:,i));
ANL_no_duct_avg(i) = mean(ANL_no_duct(:,i));
ANL_no_ice_avg(i)= mean(ANL_no_ice(:,i));

 end
 
%%
ANL_duct_noduct_avg = zeros(1,39);
ANL_duct_noice = zeros(1,39);
ANL_noduct_noice_avg = zeros(1,39);
ANL_duct_noduct = zeros(length(ANL_duct),39);
 
 
 for i=1:39
     ANL_duct_noduct_avg(i)=ANL_duct_avg(i)-ANL_no_duct_avg(:,i);
     ANL_duct_noice_avg(i)=ANL_duct_avg(i)-ANL_no_ice_avg(i);
     ANL_noduct_noice_avg(i)=ANL_no_duct_avg(i)-ANL_no_ice_avg(i);
     
     ANL_duct_noduct(:,i)=ANL_duct(:,i)-ANL_no_duct(:,i);
%      ANL_duct_noice(:,i)=ANL_duct(:,i)-ANL_no_ice(:,i);
%      ANL_noduct_noice(:,i)=ANL_no_duct_avg(:,i)-ANL_no_ice(:,i);
 end
 
 %% Figure
 
figure
plot(avg_freq,ANL_duct_noduct_avg)
hold on
plot(avg_freq,ANL_duct_noice_avg) 
plot(avg_freq,ANL_noduct_noice_avg)
hold off

ylabel(['Difference in Average ANL (dB re 1 \muPa^2/ Hz)'])
xlabel('Frequency (Hz)')
% xticks([50 100 150 200 250 300 350 400 450 500 550 600 ...
%     650 700 750 800 850 900 950 1000 1050 1100 1150 1200 1250 1300 1350 1400 1450 ...
%     1500 1550 1600 1650 1700 1750 1800 1850 1900 1950])
xticks([100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900])
xtickangle(-45)
title('Average ANL Difference at each frequency')
legend('Duct/No Duct Difference','Duct/No Ice Differnce','No Duct/No Ice Difference')

%% Figure Duct no Duct
figure
%plot(avg_freq,ANL_duct_noduct)

ylabel(['Difference in Average ANL (dB re 1 \muPa^2/ Hz)'])
xlabel('Frequency (Hz)')
% xticks([50 100 150 200 250 300 350 400 450 500 550 600 ...
%     650 700 750 800 850 900 950 1000 1050 1100 1150 1200 1250 1300 1350 1400 1450 ...
%     1500 1550 1600 1650 1700 1750 1800 1850 1900 1950])
xticks([100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900])
xtickangle(-45)
title('Average ANL Difference at each frequency')

%% Not really a point in making an imagcesc


%% covariance created not
% sigma1 = ANL_duct;
% sigma1 = ANL_no_duct;
% sigma1 = ANL_no_ice;
% Sn=zeros(size(sigma));
% Sx=zeros(size(sigma));
% for ii=1:39
%     n=v_interf_rep(:,ii)+noise_ok(:,ii);
%     s=v_interf_rep(:,ii)+noise_ok(:,ii)+v_t_rep(:,ii);
%     
%     Sn=Sn+n*n';
%     Sx=Sx+s*s';
% end
% Sn=Sn/n_snap;
% Sn_inv=inv(Sn);
% 
% Sx=Sx/n_snap;
% Sx_inv=inv(Sx);

%% Matlab generated covariance

cov_duct_freqs = cov(ANL_duct);
cov_no_duct_freqs = cov(ANL_no_duct);
cov_noice_freqs = cov(ANL_no_ice);

figure
imagesc(avg_freq,avg_freq,cov_duct_freqs)
xlabel('Frequency (Hz)')
ylabel('Frequency (Hz)')
title('Covariance between Frequencies under Ice with Duct')
c = colorbar
caxis([5 55])
axis xy

figure
imagesc(avg_freq,avg_freq,cov_no_duct_freqs)
xlabel('Frequency (Hz)')
title('Covariance between Frequencies under Ice with No Duct')
c = colorbar
caxis([5 55])
axis xy

figure
imagesc(avg_freq,avg_freq,cov_noice_freqs)
xlabel('Frequency (Hz)')
title('Covariance between Frequencies with No Ice')
c = colorbar
caxis([5 55])
axis xy

%% Correlation is NOT covariance
% correff(A) A, where the columns of A represent random variables and the rows represent observations.

% close all
%%% corr code coming from spatial
[corr_ANL_no_ice, p_no_ice] = corrcoef(ANL_no_ice); %'type','Pearson','rows','all','tail','both');
[corr_ANL_no_duct, p_no_duct] = corr(ANL_no_duct,'type','Pearson','rows','all','tail','both');
[corr_ANL_duct, p_duct] = corr(ANL_duct,'type','Pearson','rows','all','tail','both');

figure
tiledlayout(1,3)
nexttile
imagesc(avg_freq,avg_freq,corr_ANL_no_ice)
xlabel('Frequency (Hz)')
title('No Ice')
xticks([200 400 600 800 1000 1200 1400 1600 1800])
xtickangle(-45)
c = colorbar
caxis([0.4 1]);
axis xy

nexttile
imagesc(avg_freq,avg_freq,corr_ANL_no_duct)
xlabel('Frequency (Hz)')
title('Ice Without Duct')
xticks([200 400 600 800 1000 1200 1400 1600 1800])
xtickangle(-45)
c = colorbar
caxis([0.4 1]);
axis xy

nexttile
imagesc(avg_freq,avg_freq,corr_ANL_duct)
xlabel('Frequency (Hz)')
title('Ice with Duct')
xticks([200 400 600 800 1000 1200 1400 1600 1800])
xtickangle(-45)
c = colorbar
caxis([0.4 1]);
axis xy