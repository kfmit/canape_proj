%%% to run matlab from shell
%%% /usr/local/bin/matlab -nodisplay -nosplash
%%%%% and with nohup
%%% nohup /usr/local/bin/matlab -nodisplay -nosplash -r a_compute_PSD_beamform_a_la_Kinda > matlab_nohup_PSD_beamform.out &


clc
close all
clear all

format short

%%% beamforming angles in degree
%%% NB 0 is vertical/endfire (up), 90 is horizontal/broadside, 180 is
%%% vertical/endfire (down)
theta=[0 45 90 135 180];
c=1450; %%% water sound speed (m/s)
depth=[0 2.5 5 7.5];  %%% VLA depth (m)
%% user-defined parameters : LTSA option
%%% decimation factor for LTSA computation
n_decimate=1;
%%% NB: if n_decimate = 1, then no decimation 

%%% Long term spectrogram parameters
nFFT_sec= 62.5*0.001; % short window, in seconds
nIntegWin_sec= 7*60; % long window, in seconds

nOverlapFFT_sec = nFFT_sec/2; % in second
nOverlapIntegWin_sec = 0; % in seconds

step_IntegWin = 1; % One IntegWin every step_IntegWin file will be processed to compute the LTSA. 
%%%% NB: if step_IntegWin = 1, then all IntegWin will be processed

%%%% Spectrogram computation "a la Kinda"
%%%% (short term FFTs are sorted over a long term period and one keeps a percentile) 
percentil_kinda = [5 10 15 50];

nmax_IntegWin_per_DataFile = Inf;% Max number of IntegWin per data files
%%%% NB: if nmax_IntegWin_per_file = Inf, then LTSA will be computed over
%%%% entire data files.
step_DataFile = 1; % One data files every step_DataFile file will be processed to compute the LTSA. 
%%%% NB: if step_DataFile = 1, then all data files are processed

NberMaxFile=Inf; % Number of the last data file that will be processed
%%%% NB: if NberMaxFile = Inf, then all files are processed


%% Folders with acoustics data and routines

fixed_gain = 26; 
hydrophone_sensitivity=170;
chn =0:3;
Nchan=length(chn);

%%% Folder with SHRU routines
addpath('/home/julien/ju_boulot/matlab/SHRU_programs/') %% local on moliere
% addpath('/home/jbonnel/Desktop/SBCX_SHRUS/SHRU_programs/') %% local on laptop
% addpath('/home/julien/SHRU_programs/') %% on chaos
%%% Folder with SHRU data (not local, must be mounted first)
%         SHRU_Data= '/home/julien/ju_boulot/canape_chaos/Canape2016_shru_data/'; %% on chaos from moliere
%         SHRU_Data= '/home/jbonnel/canape_chaos/Canape2016_shru_data/'; %% on chaos from laptop
%         SHRU_Data= '/data3/canape/Canape2016_shru_data/'; %% on chaos from chaos
SHRU_Data= '/raid0/canape/'; %% local on moliere
%%% SHRU name/number (name of folder in SHRU_Data without suffix)
SHRU_SN='CANAPE_SHRU1_903' ;
% SHRU_SN='CANAPE_SHRU2_906' ;
% SHRU_SN='CANAPE_SHRU4_905' ;
% SHRU_SN='CANAPE_SHRU5_914' ;

SHRU_Data_dir1 = [SHRU_Data SHRU_SN '_0'];
SHRU_Data_dir2 = [SHRU_Data SHRU_SN '_1'];

%%% Read SHRU info data
%         Fileinfo = [SHRU_SN '_onchaos_frommoliere_FileInfo.mat'];
% Fileinfo = [SHRU_SN '_onchaos_fromlaptop_FileInfo.mat'];
% Fileinfo = [SHRU_SN '_onchaos_fromchaos_FileInfo.mat'];
Fileinfo = [SHRU_SN '_onmoliere_frommoliere_FileInfo.mat'];
SHRU_File = '*.D1*';  % ref SHRU folder and files files

t=0;
dt=0;
drift=0;
yr0 = 2016;

if exist(Fileinfo,'file')
    load(Fileinfo,'SHRU')
else
    %%% get the file info into the 2 data folders
    SHRU1 = sub_SHRUFileInfo(SHRU_Data_dir1, SHRU_File, t, dt,drift, yr0);  
    SHRU2 = sub_SHRUFileInfo(SHRU_Data_dir2, SHRU_File, t, dt, drift, yr0);  
    %%% SHRU1 and SHRU2 are concatenated into a single SHRU variable
    SHRU.filenames = [SHRU1.filenames; SHRU2.filenames];
    SHRU.timestamps = SHRU1.timestamps; 
    SHRU.timestamps(:,:,size(SHRU1.timestamps,3)+[1:size(SHRU2.timestamps,3)]) = SHRU2.timestamps;
    SHRU.rhfs = SHRU1.rhfs; 
    SHRU.timestamps_orig = SHRU1.timestamps_orig; 
    SHRU.timestamps_orig(:,:,size(SHRU1.timestamps_orig,3)+[1:size(SHRU2.timestamps_orig,3)]) = SHRU2.timestamps_orig;
    SHRU.rhfs_orig = SHRU1.rhfs_orig; 
    save(Fileinfo,'SHRU')
end
Nfile=size(SHRU.filenames,1);
fs_orig=SHRU.rhfs;

  

%% code
 

path_auxData = ['.' filesep 'auxData' filesep];
path_Codes = ['.' filesep 'Codes' filesep];
path_ArchivedPSDcomputation = ['.' filesep 'ArchivedPSDcomputation' filesep];

addpath(genpath(path_Codes))


vPSD_kinda=[];
vPSD_pwelch_kinda=[];
timestamp_wavDataFiles=[];
ii=1;

tic


fs=fs_orig/n_decimate;

nFFT = 2^(nextpow2(nFFT_sec*fs)); %%% in samples      
nOverlapFFT = 2^(nextpow2(nOverlapFFT_sec*fs)); %%% in samples  

nIntegWin = round(nIntegWin_sec*fs);  %%% in samples      
nOverlapIntegWin = round(nOverlapIntegWin_sec*fs);  %%% in samples      


for ww=1:step_DataFile:min(Nfile,NberMaxFile)
   
    filename=SHRU.filenames(ww,:);
    fid = fopen(filename,'rb','ieee-be');   % big endian for SHRU
    [drhs]=SHRU_getAllDRH(fid);
    fclose(fid);

    date=drhs(1).adate(1:10); 
    time=drhs(1).atime(1:8);
    date_str=[date ' ' time];
    tstartfile=datenum(date_str, 'mm/dd/yyyy HH:MM:SS'); %% [in days]      

    Nrec=length(drhs); %% number of records in the shru file

    [x,t0,header] = SHRU_getrawdata(filename,0:Nrec-1,chn); % in V, between -0.125 and 0.125    
    x=x*10^(170/20); % in uPa

    
    if n_decimate>1
        x_orig=x;
        clear x
        for cc=1:Nchan
            x(:,cc)=decimate(x_orig(:,cc),n_decimate);
        end
    end
    clear x_orig
 
    
    %%% remove constant offset (?)
    nx = size(x,1); 
    offset=mean(x);
    x=x-repmat(offset,size(x,1),1);            

    %% Input variables
    

    if nIntegWin>length(x)
        nIntegWin = length(x);
    end

    k = fix((nx-nOverlapIntegWin)/(nIntegWin-nOverlapIntegWin));

    LminusOverlap = nIntegWin-nOverlapIntegWin;
    xStart = 1:LminusOverlap:k*LminusOverlap;
    xEnd   = xStart+nIntegWin-1;
  
    
    ind_t=1;
    for theta_look=theta        
        %%% define steering vectors
        f=(0:nIntegWin-1)*fs/nIntegWin;
        k_beam=(2*pi*f/c);
        d=zeros(nIntegWin, Nchan);
        for cc=1:Nchan
            d(:,cc)=exp(1j*k_beam*depth(cc)*cosd(theta_look));
        end
        
        if ind_t == 1
            ii0=ii;
        end
        ii=ii0;
        for indIntegWin = 1:step_IntegWin:min(nmax_IntegWin_per_DataFile,k)    
            % xint has the size of nIntegWin
            xint = x(xStart(indIntegWin):xEnd(indIntegWin),:);

            ddd = datestr(addtodate(tstartfile,round(1000*xStart(indIntegWin)/fs),'millisecond'),'yyyymmddHHMMSS');
            if ind_t == 1
                timestamp_wavDataFiles = [timestamp_wavDataFiles ; {ddd}];
            end
            %%%%% beamform
            xint_f=fft(xint,nIntegWin,1);                        
            xint_f_beam=xint_f.*d;
            xint_t_beam=ifft(xint_f_beam,nIntegWin,1, 'symmetric');
            xint_mono=mean(xint_t_beam,2);

            %%%%% compute spectro
            [~,fPSD,~,vPSD_spectro] = spectrogram(xint_mono,nFFT,nOverlapFFT,nFFT,fs,'psd','onesided');            
            vPSD_kinda(ii,:,:,ind_t) =prctile(vPSD_spectro,percentil_kinda,2);  
            vPSD_pwelch_kinda(ii,:,ind_t) = mean(vPSD_spectro,2)';           
            ii=ii+1;
        end
        ind_t=ind_t+1;
        
    end
    
    z=toc;
    disp(['File ',num2str(ww),'/',num2str(min(Nfile,NberMaxFile)),' ; Elapsed time: ', num2str(z/60), ' min']);

end


elapsetipeComputeAcoustics=toc/3600;

did = dir([path_ArchivedPSDcomputation 'PSD*']);

namefile=strcat(SHRU_SN,'_beamform_sw_',num2str(nFFT_sec),'_lw_',num2str(nIntegWin_sec),'_osw_',num2str(nOverlapFFT_sec),'_olw_',num2str(nOverlapIntegWin_sec),'_theta');
for tt=1:length(theta)
   namefile=[namefile '_' num2str(theta(tt))];
end



save([path_ArchivedPSDcomputation 'PSD_' namefile '.mat'],'theta','vPSD_kinda','vPSD_pwelch_kinda','timestamp_wavDataFiles','fPSD','elapsetipeComputeAcoustics','percentil_kinda','-v7.3')   


nRawObs = length(timestamp_wavDataFiles);

nFFT_tab={ [num2str(nFFT_sec) ' / ' num2str(nFFT)] };
nIntegWin_tab={ [num2str(nIntegWin_sec) ' / ' num2str(nIntegWin)] };
nOverlapIntegWin_tab={ [num2str(nOverlapIntegWin_sec) ' / ' num2str(nOverlapIntegWin)] };
nOverlapFFT_tab={ [num2str(nOverlapFFT_sec) ' / ' num2str(nOverlapFFT)] };
n_decimate={n_decimate};
writetable(table(nFFT_tab,nIntegWin_tab,nOverlapFFT_tab,nOverlapIntegWin_tab,...
    nmax_IntegWin_per_DataFile,step_DataFile,...
    elapsetipeComputeAcoustics,nRawObs,n_decimate),[path_ArchivedPSDcomputation 'metadataAcousticComputation_' namefile '.csv'])



disp('Done')

if ~usejava('desktop')
	disp('Done')
    exit
end
