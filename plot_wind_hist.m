close all
clear all
clc


load ANL_SHRU5.mat

wind=T_ecmwf.W10;
t_wind=timestamp_num_ecmwf;





%%
clc

N_month=9;
m0=10;

hf=figure('units','normalized','outerposition',[0 0 .5 1]);
for mm=1:N_month
%     for mm=1
    t0(mm)=datenum([2016,m0+mm,1]);
    t1(mm)=datenum([2016,m0+mm,31]);
    disp(['from ', datestr(t0(mm)), ' to ' datestr(t1(mm))])
     
    ind_t=find(t_wind>=t0(mm) & t_wind<=t1(mm));
    
    %%% restric to month of interest
    wind0=wind(ind_t);
   
    subplot(3,3,mm);
    histogram(wind0,0:15)
    grid on
    xlabel('Wind speed (m/s)')
    title([datestr(t0(mm)), ' to ' datestr(t1(mm))])
    
    
%     ha(mm)=subplot(3,3,mm);
%     plot(wind_ok(ice_ok>ice_limit),anl_ok(ice_ok>ice_limit),'xr')
%     hold on
%     plot(wind_ok(ice_ok<ice_limit),anl_ok(ice_ok<ice_limit),'ob')
%     hold on
%     plot(wind_lin,anl_lin,'k', 'linewidth',2)
%     xlim([0 15])
%     ylim([60 95])
%     grid on
%     xlabel('Wind speed (m/s)')
%     ylabel('ANL (dB re 1\muPa^2)')
%     title([datestr(t0(mm)), ' to ' datestr(t1(mm))])
%     


end


%%

Nice=length(timestamp_num_ssmi);
%%%%% Period with ice
toto=find(T_ssmi.icefrac>85);
beg_ice=timestamp_num_ssmi(toto(1));
end_ice=timestamp_num_ssmi(toto(end));


ind_ice=find(t_wind>beg_ice & t_wind<end_ice);


%%%%% Period without ice
toto=find(T_ssmi.icefrac(1:floor(Nice/2))<15);
end_no_ice=timestamp_num_ssmi(toto(end));
ind_no_ice_1=find(t_wind < end_no_ice);

toto=find(T_ssmi.icefrac(end:-1:floor(Nice/2))<15);
beg_no_ice=timestamp_num_ssmi(Nice-toto(end));
ind_no_ice_2=find(t_wind > beg_no_ice);

ind_no_ice=[ind_no_ice_1 ; ind_no_ice_2];

wind_ice=wind(ind_ice);
wind_no_ice=wind(ind_no_ice);

wind_duct=wind_ice(1:end/2);
wind_noduct=wind_ice(end/2+1:end);

bins=0:15;
figure
histogram(wind_duct, bins,'Normalization', 'pdf')
hold on
histogram(wind_noduct, bins,'Normalization', 'pdf')
legend('Wind during ice and duct', 'Wind during ice wihtout duct')
xlabel('Wind speed (m/s)')
grid on

% [h, p] = ttest2(wind_duct,wind_noduct,  'Vartype','unequal')

wind_duct_m=mean(wind_duct)
wind_noduct_m=mean(wind_noduct)