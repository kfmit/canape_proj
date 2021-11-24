%%% to run matlab from shell
%%% /usr/local/bin/matlab -nodisplay -nosplash
%%%%% and with nohup
%%% nohup /usr/local/bin/matlab -nodisplay -nosplash -r b_study_soudnscape > matlab_nohup_soundcape_ice.out &

clc
% close all
clear all

format short



%% Code options
%%%% name of the LTSA to study
% namefile='CANAPE_SHRU1_903_sw_0.0625_lw_420_osw_0.03125_olw_0';
% namefile='CANAPE_SHRU2_906_sw_0.0625_lw_420_osw_0.03125_olw_0';
% namefile='CANAPE_SHRU4_953_sw_0.0625_lw_420_osw_0.03125_olw_0';
namefile='CANAPE_SHRU5_914_sw_0.0625_lw_420_osw_0.03125_olw_0';


% namedir='sw_1_lw_200_osw_0_olw_0';   %%% Roth
% namedir='sw_2_lw_300_osw_1_olw_0';   %%% Menze
% namedir='sw_60_lw_3600_osw_30_olw_0';   %%% Long
% (there must be a PSD_namedir.mat file in the ArchivedPSDcomputation
% folder)

check_for_existing_study = 0;
save_matlab_fig=0;

%% Soundscape option
plotLTSA=1;
computeEPD=1;
plotScatt=1;
computeCorr=1;
computeSpatialCorr = 0;
computePCA=0;

maxTimestampMatching = 3600 * 24; % maximum interval (in sec) to match an acoustic data with an env data


namesPlot={'W10','tp','T2062','T2061','T2075', 'Tsurf'};
% namesPlot={'W10','tp','T2085','T2073','T2050','sbe3125', 'Tsurf'};
% namesPlot={'W10','tp'};

% vectorFreqBand = [30,80,10,500,500,1000];
% vectorFreqBand = [30,80,10,500];
% vectorFreqBand = [30,80, 50,500, 500,1000, 1000,2000];
vectorFreqBand = [30,80, 250, 350, 500,1000, 1000,2000];
% vectorFreqBand = [250,350];


%%%% set up a time period for the soundscape study
time_limit=1; %%% 0 for no specific period, 1 for the period defined by t_beg and t_end
t_beg=[2017,01,15,00,00,00];  %%% date vector [YY, MM, DD, h, min, sec]
t_end=[2017,01,31,00,00,00];   %%% same

% t_beg=[2016,10,20,00,00,00];  %%% date vector [YY, MM, DD, h, min, sec]
% t_end=[2017,08,01,00,00,00];   %%% same

%%%% db range for plots (LTSA and empirical densities)
% db_lim_plot=[25 90];
db_lim_plot=[];
%%%% db_lim_plot=[] to let Matlab choose the limit based on data



%% build folder result and load acoustic features
tic
namedir=[namefile '_' datestr(t_beg,'yyyymmdd') '_' datestr(t_end,'yyyymmdd')];

name_shru=namefile(8:12);
path_auxData = ['.' filesep 'auxData_' name_shru filesep];
path_Results = ['.' filesep 'Results' filesep];
path_Codes = ['.' filesep 'Codes' filesep];
path_ArchivedPSDcomputation = ['.' filesep 'ArchivedPSDcomputation' filesep];

addpath(genpath(path_Codes))

if check_for_existing_study
    if exist([path_Results namedir],'dir') == 7
        disp('Soundscape study already done (see result folder)')
        if ~usejava('desktop')
            disp('Done')
            exit
        end
        return       
    end
end

mkdir([path_Results namedir])

copyfile([path_ArchivedPSDcomputation 'PSD_' namefile '.mat'],[path_Results namedir filesep 'PSD.mat'])

copyfile([path_ArchivedPSDcomputation 'metadataAcousticComputation_' namefile '.csv'],[path_Results namedir filesep 'metadataAcousticComputation.csv'])
copyfile(['.' filesep 'Codes'],[path_Results namedir filesep 'Codes'])

path_Results = [path_Results namedir filesep];

path_ResultsFigures = [path_Results 'Figures'];    
mkdir(path_ResultsFigures)


load([path_Results 'PSD.mat']); 


Nband=length(vectorFreqBand)/2;


timestamp_num_spectro=datenum(timestamp_wavDataFiles,'yyyymmddHHMMSS'); %%% date as a numeric array (in days since Jan 0, 0000)

%% restric acoustic data to the period of interest

if time_limit
    t_beg_num=datenum(t_beg);
    t_end_num=datenum(t_end);

    ind_t_ok=find(timestamp_num_spectro > t_beg_num & timestamp_num_spectro < t_end_num);

    timestamp_num_spectro=timestamp_num_spectro(ind_t_ok);
    vPSD_pwelch_kinda=vPSD_pwelch_kinda(ind_t_ok,:,:);
    vPSD_kinda=vPSD_kinda(ind_t_ok,:,:,:);

end



Nt_spectro=size(vPSD_pwelch_kinda,1);
Nper=length(percentil_kinda);

%% Loop over percentiles and RMS

for pp=1:Nper+1

    if pp==Nper+1
        vPSD_pp=vPSD_pwelch_kinda;
        name_pp='_rms';
    else
        vPSD_pp=squeeze(vPSD_kinda(:,:,pp,:));
        name_pp=['_p' num2str(percentil_kinda(pp))];
    end
    
%% fusion of PSD with multiple ocean observation sources        
    disp(['Merging Acoustics and Environnemental data, percentile ' ... 
        num2str(pp) '/' num2str(Nper+1)])

    if size(vPSD_pp,1)<20
        error('Not enough obs for soundscape !')
    end

    tstart=clock;

    auxDataFiles = dir([path_auxData 'variables*.csv']);

    auxData_t_psd=[]; %%% Matrix of env data interpolated on the spectro timestamp
    ecartTime_tot=[];
    auxVarNames=[];
    %%%% Nearest neighbors interpolation of the env data onto the
    %%%% spectrogram timestamp
    for aa = 1:size(auxDataFiles,1)

        if ~computeSpatialCorr && strcmp(auxDataFiles(aa).name(1:19),'variables_icemotion' )
            % don't include the icemotion variable if computeSpatialCorr =0
        else

            T = readtable([path_auxData auxDataFiles(aa).name]);
            timestamp_num_env=datenum(num2str(T.timestamp),'yyyymmddHHMMSS');

            ecartTime=zeros(Nt_spectro,1);
            vecind=zeros(Nt_spectro,1);
            for tt=1:Nt_spectro
                [ecartTime(tt),vecind(tt)] = min(abs(timestamp_num_env - timestamp_num_spectro(tt))); %% in days
            end   
            ecartTime=ecartTime/3600; %% in sec

            auxData_t_psd = [auxData_t_psd , T{vecind,2:end}];
            ecartTime_tot = [ecartTime_tot , repmat(ecartTime,1,size(T{vecind,2:end},2))];
            auxVarNames=[auxVarNames , T.Properties.VariableNames(2:end)];

        end
    end


    % put NaN for time observation that do not fit within maxTimestampMatching      
    auxData_t_psd(ecartTime_tot> maxTimestampMatching)=NaN; 

    % remove auxiliary variables with too many NaN values
    removeVar=find(sum(isnan(auxData_t_psd),1) > 0.75*size(auxData_t_psd,1));
    nRemovedAuxVar = length(removeVar);
    auxData_t_psd(:,removeVar)=[];
    auxVarNames(removeVar)=[];

    elapsetipeMergingAux = etime(clock,tstart);

    disp(['Acoustic and Environmental data merged in ', num2str(elapsetipeMergingAux), ' sec ; percentile ' ... 
        num2str(pp) '/' num2str(Nper+1)])

    cpt=1;
    for lll=1:length(namesPlot)
        ii=find(~cellfun(@isempty,(strfind(auxVarNames,namesPlot{lll}))),1,'first');
        if ~isempty(ii), vecAuxVarToPlot(cpt)=ii; cpt=cpt+1; end
    end




    %% feature matrix    
    Nchan=length(chn);

    for cc=1:Nchan
        vPSD=vPSD_pp(:,:,cc);


        vecSPL=zeros(Nt_spectro,Nband);

        cpt=1;
        for bb=1:2:2*Nband
           f1=vectorFreqBand(bb);
           f2=vectorFreqBand(bb+1);
           vecSPL(:,cpt)= 10*log10(mean(vPSD(:,fPSD>f1 & fPSD<f2),2));  
           nameSPL{cpt}=strcat('SPL\_',num2str(f1),'\_',num2str(f2));
           cpt=cpt+1;
        end


        %%% on concatene les SPL et les variables environnementales
        FeatureMatrix = [vecSPL,auxData_t_psd(:,vecAuxVarToPlot)];

        %%% on garde uniquement les indices temporels ou ni les SPL ni les
        %%% variables environnmentales sont des NaN
        keepingIndObs = unique(find(sum(isnan(FeatureMatrix),2) < 1));      
        nMatchedStandardized = length(keepingIndObs);
        FeatureMatrix_noNan = FeatureMatrix(keepingIndObs,:);

        %%% on normalise chaque variable (moyenne nulle et variance unitaire)
        FeatureMatrix_noNan_norma = zscore(FeatureMatrix_noNan,[],1);

        names_statSoundscapeFeatureMatrix = [nameSPL, auxVarNames(vecAuxVarToPlot)];  



    %% long-term spectro + aux data
        if plotLTSA
    
        disp(['Plotting LTSA, channel ', num2str(cc), '/', num2str(Nchan), ...
            ', percentile ', num2str(pp), '/', num2str(Nper+1)])
        Opt.OffsetFreqDescriptors=40;
        Opt.InterFreqDescriptors=200;
        Opt.NberLabelX=20;
        Opt.TimeStampFormat='yyyy/mm/dd:HH';
        Opt.MV_Apply_MedFilt = 40;
        Opt.AuxVariableNames = auxVarNames(vecAuxVarToPlot);

        figure('visible','off');
        Plot_LongTermAverageSpectro_jb(10*log10(vPSD),fPSD,timestamp_num_spectro,auxData_t_psd(:,vecAuxVarToPlot),Opt, db_lim_plot)
        FormatFigures([path_ResultsFigures filesep 'LTSA_chn_' num2str(chn(cc)) name_pp], save_matlab_fig)

        disp(['LTSA saved, channel ', num2str(cc), '/', num2str(Nchan), ...
            ', percentile ', num2str(pp), '/', num2str(Nper+1)])

        end
    %% empirical proba density   
        if computeEPD
        
        disp(['Computing spectrum probability density, channel ', num2str(cc), '/', num2str(Nchan), ...
            ', percentile ', num2str(pp), '/', num2str(Nper+1)])
        figure('visible','off');
        A=10*log10(vPSD);        
        prctilePlot_jb;
        FormatFigures([path_ResultsFigures filesep 'EPD_chn_' num2str(chn(cc)) name_pp], save_matlab_fig)
        disp(['Spectrum probability density saved, channel ', num2str(cc), '/', num2str(Nchan), ...
            ', percentile ', num2str(pp), '/', num2str(Nper+1)])

        end
    %% scattering plots  
       
         Vec_IndexResponse = [1:Nband];
         Vec_IndexExplain = [Nband+1: Nband + length(vecAuxVarToPlot)];
    
         if plotScatt
        
        disp(['Making scatter plots, channel ', num2str(cc), '/', num2str(Nchan), ...
            ', percentile ', num2str(pp), '/', num2str(Nper+1)])

        LabelLeg='off';

        figure('visible','off');
        gplotmatrix_jb(FeatureMatrix,Vec_IndexResponse,names_statSoundscapeFeatureMatrix);
        FormatFigures([path_ResultsFigures filesep 'scatAcoust_chn_' num2str(chn(cc)) name_pp], save_matlab_fig)

        if pp==1 && cc==1
            figure('visible','off');
            gplotmatrix_jb(FeatureMatrix,Vec_IndexExplain,names_statSoundscapeFeatureMatrix);
            FormatFigures([path_ResultsFigures filesep 'scatEnv'], save_matlab_fig)
        end

        disp(['Scatter plots saved, channel ', num2str(cc), '/', num2str(Nchan), ...
            ', percentile ', num2str(pp), '/', num2str(Nper+1)])

        end
    %% scattering and pearson correlation coefs
        if computeCorr
    
        disp(['Computing correlations, channel ', num2str(cc), '/', num2str(Nchan), ...
            ', percentile ', num2str(pp), '/', num2str(Nper+1)])
        %%% If removeOutliers=1 then outliers will be removed from correlation
        %%% analysis
        removeOutliers=0;

        Matrix_X=FeatureMatrix_noNan_norma(:,Vec_IndexExplain);
        Matrix_Y=FeatureMatrix_noNan_norma(:,Vec_IndexResponse);   

        %parameters for figure and panel size
        subplotsx=length(Vec_IndexExplain);
        subplotsy=length(Vec_IndexResponse);  

        figure('visible','off');   
        SubplotCorr;       
        FormatFigures([path_ResultsFigures filesep 'PairWiseCorr_chn_' num2str(chn(cc)) name_pp], save_matlab_fig)

        disp(['Correlations saved, channel ', num2str(cc), '/', num2str(Nchan), ...
            ', percentile ', num2str(pp), '/', num2str(Nper+1)])

        end
    %% spatial Pearson correlation with ice motion  

        cpt=1;
        for bb=1:2:2*Nband
           f1=vectorFreqBand(bb);
           f2=vectorFreqBand(bb+1);
           vecSPL(:,cpt)= 10*log10(mean(vPSD(:,fPSD>f1 & fPSD<f2),2));  
           nameSPL{cpt}=strcat('SPL\_',num2str(f1),'\_',num2str(f2));
           cpt=cpt+1;
        end


        if computeSpatialCorr
%                 T=readtable([path_auxData 'info_variables_icemotion.csv']);
            T=readtable([path_auxData 'info_variables_icemotion_osisaf.csv']);
            
            for ff=1:Nband
                f1=vectorFreqBand(2*(ff-1)+1);
                f2=vectorFreqBand(2*(ff-1)+2);

                disp(['Computing spatial correlations, channel ', num2str(cc), '/', num2str(Nchan), ...
                    ', freq ' num2str(f1), '-', num2str(f2), ' Hz', ', percentile ', num2str(pp), '/', num2str(Nper+1)])
                indSpatialCorr=find(~cellfun(@isempty,(strfind(auxVarNames,'col'))));

                vecR=NaN(length(indSpatialCorr),3);
                vecP=NaN(length(indSpatialCorr),1);

                for ii=1:length(indSpatialCorr)
                                                           
                    X = [vecSPL(~isnan(auxData_t_psd(:,indSpatialCorr(ii))),ff), auxData_t_psd(~isnan(auxData_t_psd(:,indSpatialCorr(ii))),indSpatialCorr(ii))];
                    X_norma=zscore(X,[],1);
                    [R,PValue] = corr(X_norma(:,1),X_norma(:,2),'type','Pearson','rows','all','tail','both');

                    indCol=find(~cellfun(@isempty,(strfind(T.col_number,auxVarNames{indSpatialCorr(ii)}))),1,'first');        
                    vecR(ii,:)=[R,T.ind_lat(indCol),T.ind_lon(indCol)];
                    vecP(ii)=PValue;       
                end

%                 plot_spatial_corr;
                plot_spatial_corr_osisaf;
                FormatFigures([path_ResultsFigures filesep 'SpatialPearson_chn_'  num2str(chn(cc)) ...
                    '_freq_' num2str(f1) '-' num2str(f2) name_pp], save_matlab_fig)

                disp(['Spatial correlations saved, channel ', num2str(cc), '/', num2str(Nchan), ...
                    ', freq ' num2str(f1), '-', num2str(f2), ' Hz', ', percentile ', num2str(pp), '/', num2str(Nper+1)])

            end

        end


    %% PCA   

        if computePCA
            disp(['Computing PCA, channel ', num2str(cc), '/', num2str(Nchan), ...
                ', percentile ', num2str(pp), '/', num2str(Nper+1)])
            [wcoeff,score,latent,tsquared,explained] = pca(FeatureMatrix_noNan_norma,'VariableWeights','variance');
            %Transform the coefficients so that they are orthonormal
            coefforth = diag(std(FeatureMatrix_noNan_norma))\wcoeff;

            figure('visible','off');
            biplot(coefforth(:,1:2),'scores',score(:,1:2),'varlabels',names_statSoundscapeFeatureMatrix);
            grid on 
            FormatFigures([path_ResultsFigures filesep 'PCA_chn_' num2str(chn(cc)) name_pp], save_matlab_fig)

            disp(['PCA saved, channel ', num2str(cc), '/', num2str(Nchan), ...
                ', percentile ', num2str(pp), '/', num2str(Nper+1)])
        end   
    end
end



%% End
elapsetipeSoundscapeWorkflow = toc / 60;
maxTimestampMatching_hour=maxTimestampMatching / 3600 ;
writetable(table(names_statSoundscapeFeatureMatrix,nMatchedStandardized,maxTimestampMatching_hour,elapsetipeMergingAux,elapsetipeSoundscapeWorkflow),[path_Results 'metadataSoundscapeWorkflow.csv'])



disp('Done')

if ~usejava('desktop')
	disp('Done')
    exit
end
