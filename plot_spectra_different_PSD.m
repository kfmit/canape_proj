close all
clear all
clc


color={'r', 'g', 'k', 'b'};

ind_p=3;
ind_capt=1;

plot_EDP=0;


for ii=1:4
    if ii==4
        load ArchivedPSDcomputation/PSD_CANAPE_SHRU5_914_sw_0.0625_lw_420_osw_0.03125_olw_0.mat
        ind_cut=34;   
        param_ltsa='sw 0.0625 lw 420';
    elseif ii==3
        load ArchivedPSDcomputation/PSD_CANAPE_SHRU5_914_sw_1_lw_200_osw_0_olw_0.mat
        ind_cut=1; 
        param_ltsa='sw 1 lw 200';
    elseif ii==2
        load ArchivedPSDcomputation/PSD_CANAPE_SHRU5_914_sw_2_lw_300_osw_1_olw_0.mat
        ind_cut=1;
        param_ltsa='sw 2 lw 300';
    elseif ii==1
        load ArchivedPSDcomputation/PSD_CANAPE_SHRU5_914_sw_60_lw_3600_osw_30_olw_0.mat
        ind_cut=2;
        param_ltsa='sw 60 lw 3600';
    end
    
    vPSD_kinda(1:ind_cut,:,:,:)=[];
    vPSD_pwelch_kinda(1:ind_cut,:,:)=[];
    timestamp_wavDataFiles(1:ind_cut)=[];

    T = readtable('auxData_SHRU5/variables_ASI-SSMI_SHRU5.csv');
    timestamp_num_env=datenum(num2str(T.timestamp),'yyyymmddHHMMSS');

    timestamp_num_spectro=datenum(timestamp_wavDataFiles,'yyyymmddHHMMSS');


    ANL=squeeze(vPSD_kinda(:,:,ind_p,ind_capt));

    LTSA=squeeze(vPSD_pwelch_kinda(:,:,ind_capt));
    
    %% Define ice / no-ice

    p_no_ice=15;
    ind_no_ice = find(T.icefrac<p_no_ice);
    [~, toto]=max(abs(diff(ind_no_ice)));

    p_ice=85;
    ind_ice=find(T.icefrac>p_ice);

    beg_ice=timestamp_num_env(ind_ice(1));
    end_ice=timestamp_num_env(ind_ice(end));
    ind_acoust_ice=find(timestamp_num_spectro>beg_ice & timestamp_num_spectro<end_ice);

    end_no_ice_acoust=find(timestamp_num_spectro<timestamp_num_env(ind_no_ice(toto-1)));
    end_no_ice_acoust=end_no_ice_acoust(end);
    beg_no_ice_acoust=find(timestamp_num_spectro>timestamp_num_env(ind_no_ice(toto+1)));
    beg_no_ice_acoust=beg_no_ice_acoust(1);
    ind_acoust_no_ice=[1:end_no_ice_acoust beg_no_ice_acoust:length(timestamp_num_spectro)];


    figure, imagesc(timestamp_num_spectro, fPSD, 10*log10(LTSA).'), axis xy, datetick2('x')
    figure, imagesc(timestamp_num_spectro(ind_acoust_no_ice), fPSD, 10*log10(LTSA(ind_acoust_no_ice,:)).'), axis xy


    %% Plot spectra

    p_plot=[1 10 25 50 75 90 99];

    perc_ANL_no_ice=prctile(ANL(ind_acoust_no_ice,:),p_plot, 1);
    perc_ANL_ice=prctile(ANL(ind_acoust_ice,:),p_plot, 1);
    perc_LTSA_no_ice=prctile(LTSA(ind_acoust_no_ice,:),p_plot, 1);
    perc_LTSA_ice=prctile(LTSA(ind_acoust_ice,:),p_plot, 1);


    figure(1)
    p1=subplot(221);
    semilogx(fPSD, 10*log10(perc_ANL_no_ice), 'Color', color{ii})
    grid on
    title('No ice - ANL')
%     legend('1 %', '10 %', '25 %', '50 %', '75 %', '90 %', '99 %')
    xlabel('Frequency (Hz)')
    ylabel('dB / Hz ref 1 \muPa')
    hold on

    p2=subplot(222);
    semilogx(fPSD, 10*log10(perc_ANL_ice), 'Color', color{ii})
    grid on
    title('Ice covered - ANL')
%     legend('1 %', '10 %', '25 %', '50 %', '75 %', '90 %', '99 %')
    xlabel('Frequency (Hz)')
    ylabel('dB / Hz ref 1 \muPa')
    hold on

    p3=subplot(223);
    semilogx(fPSD, 10*log10(perc_LTSA_no_ice), 'Color', color{ii})
    grid on
    title('No ice - LTSA')
%     legend('1 %', '10 %', '25 %', '50 %', '75 %', '90 %', '99 %')
    xlabel('Frequency (Hz)')
    ylabel('dB / Hz ref 1 \muPa')
    hold on

    p4=subplot(224);
    semilogx(fPSD, 10*log10(perc_LTSA_ice), 'Color', color{ii})
    grid on
    title('Ice covered - LTSA')
%     legend('1 %', '10 %', '25 %', '50 %', '75 %', '90 %', '99 %')
    xlabel('Frequency (Hz)')
    ylabel('dB / Hz ref 1 \muPa')
    hold on

    linkaxes([p1,p2,p3,p4],'xy')
    
    %%%% plot EDP
    if plot_EDP
        hind = 2;         
        dbvec = 20:hind:120;
        A=10*log10(LTSA(ind_acoust_no_ice,:));

        %%% Compute probability densities
        d = hist(A,dbvec)/(hind*size(A,1));  %SPD array
        d(d == 0) = NaN;   %suppress plotting of empty hist bins

        figure(2)
        subplot(2,2,ii)
        g = pcolor(fPSD,dbvec,d); 
        set(g,'LineStyle','none')
        % set(gca, 'XScale', 'log')
        colorbar
        grid on
        title(param_ltsa)
        caxis([0.02 0.2])
        
    end

    
    
end
figure(1)
xlim([1 2000])
ylim([20 120])



%%