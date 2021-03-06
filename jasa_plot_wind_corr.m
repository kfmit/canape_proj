close all
clear all
clc


load ANL_SHRU5_newfreq.mat
% load ANL_SHRU5.mat

wind=T_ecmwf.W10;
t_wind=timestamp_num_ecmwf;
t_anl=timestamp_num_spectro;

% f1=[40 450 900 1250 250];
% f2=[60 550 1100 1750 350];
ff = 4;
anl=SPL_ANL(:,ff);   %%% focus on the band 250-350 Hz, column 5

ice_frac=T_ssmi.icefrac;
t_ice=timestamp_num_ssmi;
ice_limit=50;

ice_frac_interp=interp1(t_ice,ice_frac,t_wind);


%% average the anl
Nwind=length(wind);
anl_ave=zeros(size(wind));
anl_ave2=zeros(size(wind));


for tt=1:Nwind-1
    [~, b1]=min(abs(t_anl-t_wind(tt)));
    [~, b2]=min(abs(t_anl-t_wind(tt+1)));
    db=b2-b1;
    if db>10
        i_min=b1-round(db/2);
        i_min2=b1-round(db);
        if i_min<1
            i_min=1;
        end
        if i_min2<1
            i_min2=1;
        end
        i_max=min([length(t_anl),b1+round(db/2)]);
        i_max2=min([length(t_anl),b1+round(db)]);
        anl_ave(tt)=mean(anl(i_min:i_max));
        anl_ave2(tt)=mean(anl(i_min2:i_max2));
    end
end

% figure,
% plot(t_anl,anl)
% hold on
% plot(t_wind,anl_ave)
% hold on
% plot(t_wind,anl_ave2)




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
    anl0=anl_ave(ind_t);
    t0=t_wind(ind_t);
    ice0=ice_frac_interp(ind_t);
    
    %%% keep only good points
    ind_keep=find(anl0>20);
    wind_ok=wind0(ind_keep);
    anl_ok=anl0(ind_keep);
    t_ok=t0(ind_keep);
    ice_ok=ice0(ind_keep);
    
    tbl = table(wind_ok,anl_ok,'VariableNames', {'W10','ANL'});
    
    mdl = fitlm(tbl,'ANL ~ W10');
    R(mm)=sqrt(mdl.Rsquared.Ordinary);
    Pval(mm)=mdl.Coefficients.pValue(2);
    
    wind_lin=linspace(0,20,10);
    anl_lin=mdl.Coefficients.Estimate(1)+mdl.Coefficients.Estimate(2)*wind_lin;
      
    
    ha(mm)=subplot(3,3,mm);
    plot(wind_ok(ice_ok>ice_limit),anl_ok(ice_ok>ice_limit),'xr')
    hold on
    plot(wind_ok(ice_ok<ice_limit),anl_ok(ice_ok<ice_limit),'ob')
    hold on
    plot(wind_lin,anl_lin,'k', 'linewidth',2)
    xlim([0 15])        % change as needed for each iteration
    ylim([65 102])
    grid on
    xlabel('Wind speed (m/s)')
    ylabel('ANL (dB re 1\muPa^2)')
    title({[datestr(t0(mm))], [' to ' datestr(t1(mm))]})
    
    if mm==1
        legend('IC>50%', 'IC<50%','location','northwest')
    end

end
sgtitle(['Wind Correlation for [' num2str(f1(ff)) ' - ' num2str(f2(ff)) '] Hz'])



arrayfun(@(x) pbaspect(x, [1 1 1]), ha);
drawnow;
pos = arrayfun(@plotboxpos, ha, 'uni', 0);
dim = cellfun(@(x) x.*[1 1 0.5 0.5], pos, 'uni',0);

for mm = 1:N_month
    annotation(hf, 'textbox',  dim{mm}, 'String', ['R=', num2str(round(R(mm)*100)/100)], 'vert', ...
        'bottom', 'FitBoxToText','on', 'Color', 'k', 'fontsize',17,'fontweight',...
        'bold','LineStyle','none');
end

%% limits for each band

% for 250-350
% xlim([0 15])        % change as needed for each iteration
% ylim([60 95])

% for 40-60
% xlim({})
% ylim([60 95])

% for 459-550
% xlim({})
% ylim([60 95])

% for 900-1100
% xlim({})
% ylim([60 95])

% for 1250-1750
