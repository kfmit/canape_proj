close all
clear all
clc

%% Temperature
load temp_shru5.mat

figure
plot(tm0,temp_shru_ok)
datetick('x','keeplimits');
xlabel('Date');
ylabel('Temperature (deg C)')
grid on
legend('z=160 m', 'z=300 m', 'z=445 m')

temp_shru_norm=zeros(size(temp_shru_ok));
for nn=1:3
    toto=temp_shru_ok(:,nn);
    mu=nanmean(toto);
    sigma=nanstd(toto);
    temp_shru_norm(:,nn)=(toto-mu)/sigma;
end


%% Acoustics
load ANL_SHRU5.mat
t=timestamp_num_spectro;
SPL_ANL_nobeamform=SPL_ANL;

load ANL_SHRU5_beamform.mat


%%%% Limit to a single frequency band 250-350 Hz and 500-1000
ff=5; %%% 250-350 Hz
SPL_ANL_ok=SPL_ANL_nobeamform(:,ff);
ff=4; %%% 500-1000 Hz
SPL_ANL_ok2=SPL_ANL_nobeamform(:,ff);
ff=2;  %%% 250-350 Hz
SPL_ANL_up_ok=squeeze(SPL_ANL(:,ff,1));
SPL_ANL_hor_ok=squeeze(SPL_ANL(:,ff,3));
ff=5; %%% 500-1000 Hz
SPL_ANL_up_ok2=squeeze(SPL_ANL(:,ff,1));
SPL_ANL_hor_ok2=squeeze(SPL_ANL(:,ff,3));

SPL=zeros(size(SPL_ANL_ok,1),6);
SPL(:,1)=SPL_ANL_ok;
SPL(:,2)=SPL_ANL_hor_ok;
SPL(:,3)=SPL_ANL_up_ok;
SPL(:,4)=SPL_ANL_ok2;
SPL(:,5)=SPL_ANL_up_ok2;
SPL(:,6)=SPL_ANL_hor_ok2;

SPL_norm=zscore(SPL, [], 1);

figure
plot(t, SPL_norm(:,1:3))

%%
figure
for nn=1:3
    ax{nn}=subplot(3,1,nn);
    plot(tm0, temp_shru_norm(:,nn))
    grid on
    hold on
    plot(t, SPL_norm(:,1))
    datetick
    title(['Depth = ' num2str(z(nn)) ' m'])
    legend('Normalized temperature', 'Normalized SPL')
end

linkaxes([ax{1}, ax{2}, ax{3}], 'x')

%%
duct=temp_shru_ok(:,1)./temp_shru_ok(:,3);
mu=nanmean(duct);
sigma=nanstd(duct);
duct_norm=(duct-mu)/sigma;

figure,
plot(tm0, duct_norm)
hold on
plot(t, SPL_norm(:,1))
datetick
legend('Duct proxy', 'Normalized SPL')