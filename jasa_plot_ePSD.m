close all
clear all
clc



T = readtable('auxData_SHRU5/variables_ASI-SSMI_SHRU5.csv');
timestamp_num_env=datenum(num2str(T.timestamp),'yyyymmddHHMMSS');

% load ArchivedPSDcomputation/PSD_CANAPE_SHRU5_914_sw_0.0625_lw_420_osw_0.03125_olw_0.mat
load ArchivedPSDcomputation/PSD_CANAPE_SHRU5_914_sw_1_lw_200_osw_0_olw_0.mat
ind_cut=34;  
ind_capt=4;



vPSD_kinda(1:ind_cut,:,:,:)=[];
vPSD_pwelch_kinda(1:ind_cut,:,:)=[];
timestamp_wavDataFiles(1:ind_cut)=[];


timestamp_num_spectro=datenum(timestamp_wavDataFiles,'yyyymmddHHMMSS');

ind_p=3;
LTSA=squeeze(vPSD_kinda(:,:,ind_p,ind_capt));
% LTSA=squeeze(vPSD_pwelch_kinda(:,:,ind_capt));

%%
clc

N_month=9;
m0=10;

hind = 1;         
dbvec = 30:hind:90;

% originals
% min_f=25;  %% [Hz]
% max_f=350;  %% [Hz]

min_f=25;  %% [Hz]
max_f=350;  %% [Hz]

ind_f=find(fPSD<max_f & fPSD>min_f);

ePSD=zeros(N_month,length(dbvec),length(fPSD(ind_f)));
        
for mm=1:N_month
    t0(mm)=datenum([2016,m0+mm,1]);
    t1(mm)=datenum([2016,m0+mm,31]);
    disp(['from ', datestr(t0(mm)), ' to ' datestr(t1(mm))])
    
    ind_t=find(timestamp_num_spectro>=t0(mm) & timestamp_num_spectro<=t1(mm));
    
    A=10*log10(LTSA(ind_t,ind_f));
    
     %%% Compute probability densities
    d = hist(A,dbvec)/(hind*size(A,1));  %SPD array
    d(d == 0) = NaN;   %suppress plotting of empty hist bins
    
    ePSD(mm,:,:)=d;
    
    ind_t_ice=find(timestamp_num_env>=t0(mm) & timestamp_num_env<=t1(mm));
    ic(mm)=nanmean(T.icefrac(ind_t_ice));
    
end




%%
close all
clc

hf=figure('units','normalized','outerposition',[0 0 .5 1]);
for mm=1:N_month
    ha(mm)=subplot(3,3,mm);
    g = pcolor(fPSD(ind_f),dbvec,squeeze(ePSD(mm,:,:))); 
    set(g,'LineStyle','none')
    grid on
    caxis([0.02 0.2])
    ylim([dbvec(1)-5 dbvec(end)-10])
    title([datestr(t0(mm)), ' to ' datestr(t1(mm))])
    xlabel('Frequency (Hz)')
    ylabel('ANL (dB re 1 \muPa^2/Hz)')
end

c=colorbar(ha(end));
c.Location = 'southoutside'; 
c.Label.String = 'Probability density';
set(c, {'Position'}, mat2cell(vertcat(c.Position) .* [1 1 .92, 1], ones(size(c(:))),4))

arrayfun(@(x) pbaspect(x, [1 1 1]), ha);
drawnow;
pos = arrayfun(@plotboxpos, ha, 'uni', 0);
dim = cellfun(@(x) x.*[1 1 0.5 0.5], pos, 'uni',0);
for mm = 1:N_month
    annotation(hf, 'textbox',  dim{mm}, 'String', ['IC=' num2str(round(ic(mm))) '%'],...
        'vert', 'bottom', 'FitBoxToText','on', 'Color', 'r', 'fontsize',18,'fontweight','bold',...
        'LineStyle','none');
end




