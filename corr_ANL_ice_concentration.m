clear all
close all
clc

clear all
close all
clc

load ANL_SHRU5.mat
t=timestamp_num_spectro;
SPL_ANL_nobeamform=SPL_ANL;

load ANL_SHRU5_beamform.mat

t_ice=timestamp_num_ssmi;
ice=T_ssmi.icefrac;

%% Limit to a single frequency band 250-350 Hz
ff=5;
SPL_ANL_ok=SPL_ANL_nobeamform(:,ff);
ff=2;
SPL_ANL_up_ok=squeeze(SPL_ANL(:,ff,1));
SPL_ANL_hor_ok=squeeze(SPL_ANL(:,ff,3));

clear SPL_ANL SPL_ANL_raw SPL_ANL_nobeamform

figure
yyaxis left
plot(t,SPL_ANL_ok)
ylabel('ANL (dB)')
yyaxis right
plot(t_ice,ice, 'linewidth',2)
ylabel('Ice concetration (%)')
grid on
datetick('x')


%% Average acoustic on ice time_scale

t_deb=t(1);
t_end=t(end);

toto=find(t_ice>t_deb & t_ice<t_end);

ice_ok=ice(toto);
t_ice_ok=t_ice(toto);

N_ice=length(ice_ok);
SPL_ANL_ave=zeros(size(t_ice_ok));

for ii=2:N_ice-1
    t1=t_ice_ok(ii-1);
    t2=t_ice_ok(ii);
    t3=t_ice_ok(ii+1);
    
    t_deb=(t1+t2)/2;
    t_end=(t2+t3)/2;
    toto=find(t>=t_deb & t<=t_end);
    
    SPL_ANL_ave(ii)=mean(SPL_ANL_ok(toto));    
end

SPL_ANL_ave(1)=NaN;
SPL_ANL_ave(end)=NaN;

% figure
% yyaxis left
% plot(t_ice_ok,SPL_ANL_ave)
% yyaxis right
% plot(t_ice_ok,ice_ok)
% grid on
% datetick('x')

%%
x1=SPL_ANL_ave(2:end-1);
x1_norm=(x1-mean(x1))/std(x1);
x2=ice_ok(2:end-1);
x2_norm=(x2-mean(x2))/std(x2);

[rho,pval] = corr(x1,x2,'type','Pearson','rows','all','tail','both')



X = [x1_norm,x2_norm];
tbl = table(X(:,1),X(:,2),'VariableNames', {'ANL','Ice_concentration'});
mdl = fitlm(tbl,'ANL ~ Ice_concentration');
% p_plot=plotAdded(mdl);


% figure
% yyaxis left
% plot(t_ice_ok,SPL_ANL_ave)
% yyaxis right
% plot(t_ice_ok,ice_ok)
% grid on
% datetick('x')

%%
figure
yyaxis left
plot(t_ice_ok,SPL_ANL_ave, 'linewidth',2)
ylabel('ANL (dB)')
yyaxis right
plot(t_ice_ok,ice_ok, 'linewidth',2)
ylabel('Ice concetration (%)')
grid on
datetick('x')
