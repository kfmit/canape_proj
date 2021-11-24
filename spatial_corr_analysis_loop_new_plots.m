clear all
close all
clc

load ANL_SHRU5.mat
sitename = 'SHRU5';
gps_site = [72+54.4580/60 , -(157+29.2442/60)];  

load_matlab_ice = 0;


%% load env_data and interpolate on ANL time
t=timestamp_num_spectro;
if ~load_matlab_ice

    T = readtable('./auxData_SHRU5/variables_icemotion_osisaf.csv');
    timestamp_num_env=datenum(num2str(T.timestamp),'yyyymmddHHMMSS');

%     maxTimestampMatching = 3600 * 24; % maximum interval (in sec) to match an acoustic data with an env data
% 
%     ecartTime=zeros(Nt,1);
%     vecind=zeros(Nt,1);
%     for tt=1:Nt
%         [ecartTime(tt),vecind(tt)] = min(abs(timestamp_num_env - t(tt))); %% in days
%     end   
%     ecartTime=ecartTime/3600; %% in sec
%     auxData_t_psd =T{vecind,2:end};
%     % put NaN for time observation that do not fit within maxTimestampMatching      
%     auxData_t_psd(ecartTime> maxTimestampMatching)=NaN;     
auxVarNames=T.Properties.VariableNames(2:end);

    %%%%% new interp
    auxData=T{:,2:end};
    auxData_t_psd=interp1(timestamp_num_env, auxData, t);
       


     % remove auxiliary variables with too many NaN values
    removeVar=find(sum(isnan(auxData_t_psd),1) > 0.75*size(auxData_t_psd,1));
    nRemovedAuxVar = length(removeVar);
    auxData_t_psd(:,removeVar)=[];
    auxVarNames(removeVar)=[];

    T_info=readtable('./auxData_SHRU5/info_variables_icemotion_osisaf.csv');

    save auxData_icedrift_on_psd_time_new_interp auxData_t_psd  auxVarNames T_info
else
    load auxData_icedrift_on_psd_time_new_interp.mat
end
Nvar=length(auxVarNames);


load ./Codes/lat_lon_osisaf
% figure
% imagesc(1:Nvar, t,auxData_t_psd)
% datetick2('y')


%% Must loop from Nov 1st to May 31

% t_beg=[2016,12,01,00,00,00];  %%% date vector [YY, MM, DD, h, min, sec]
% t_end=[2017,02,01,00,00,00];   %%% date vector [YY, MM, DD, h, min, sec]

t0=[2016,11,01,00,00,00];
t0_num=datenum(t0);

t1=[2017,05,31,00,00,00];
t1_num=datenum(t1);

Nloop=length(t0_num:15:t1_num);
date_loop=zeros(2,Nloop);

corr_spa=NaN([size(latitude) Nf Nloop]);

tt=1;
for t_num_loop=t0_num:15:t1_num
    t_beg_num=t_num_loop;
    t_end_num=t_num_loop+60;

    disp(datestr(t_beg_num, 'yyyymmdd'))
    %% Restric to time of interest


    ind_t_ok=find(t > t_beg_num & t < t_end_num);

    t_ok=t(ind_t_ok);
    Nt=length(t);
    auxData_t_psd_ok=auxData_t_psd(ind_t_ok,:);
    SPL_ANL_ok=SPL_ANL(ind_t_ok,:);
%     SPL_raw_ok=SPL_raw(ind_t_ok,:);

    % figure
    % imagesc(1:Nvar, t_ok,auxData_t_psd_ok)
    % shading flat
    % datetick2('y')

    %% Compute correlation

    vecR=NaN(Nvar,3, Nf);
    vecP=NaN(Nvar,1, Nf);

    vecSPL=SPL_ANL_ok;

    for ff=1:Nf
        for ii=1:Nvar
            ind_no_nan=~isnan(auxData_t_psd_ok(:,ii));
            %%%% remove corr data which are nearly constant and/or nearly all
            %%%% NaN
            if length(unique(auxData_t_psd_ok(ind_no_nan,ii)))>10
                X = [vecSPL(ind_no_nan,ff), auxData_t_psd_ok(ind_no_nan,ii)];
                X_norma=zscore(X,[],1);
                [R,PValue] = corr(X_norma(:,1),X_norma(:,2),'type','Pearson','rows','all','tail','both');

                indCol=find(~cellfun(@isempty,(strfind(T_info.col_number,auxVarNames{ii}))),1,'first');        
                vecR(ii,:,ff)=[R,T_info.ind_lat(indCol),T_info.ind_lon(indCol)];
                vecP(ii,ff)=PValue;       

            end
        end
    end

    %% Plot results

    %%%% prepare plot
    latlimit=[65 85];
    dlon=40;
    lonlimit=[gps_site(2)-dlon gps_site(2)+dlon];
    centralmeridian=-160;
    parallel=[70 75 80];
    MLabelLocation=[-180:20:180];
    Mpos=latlimit(1)+2;
    PLabelMeridian=30;
    land = shaperead('landareas.shp', 'UseGeoCoords', true);
    
    
    for ff=1:Nf
        figure('visible','off');

        for ii=1:length(vecR)
            if vecP(ii,ff)<0.05
                corr_spa(vecR(ii,2,ff),vecR(ii,3,ff),ff,tt)=abs(vecR(ii,1,ff)); 
            elseif vecP(ii,ff)>=0.05
                corr_spa(vecR(ii,2,ff),vecR(ii,3,ff),ff,tt)=NaN; 
            end
        end

        %%% set axes
        maph=axesm('MapProjection','lambertstd','MapLatLimit',latlimit,'MapLonLimit',lonlimit);

        %%% add grid
        setm(maph,'Grid','on','Glinewidth',1,'PLineLocation',parallel, 'MLineLocation',MLabelLocation);

        %%% add grid labeling
        setm(maph,'Fontangle','normal',...
          'FontSize',12,'fontweight','b',...
          'MeridianLabel','on',...
          'MLabelLocation',MLabelLocation,...
          'MLabelParallel',Mpos,...
          'ParallelLabel','on',...
          'PLabelLocation',parallel,...
          'PLabelMeridian',PLabelMeridian);

        %%% add land
        geoshow(maph, land, 'FaceColor',[0.80 0.80 0.80],'EdgeColor',0.30*[1 1 1]);

        %%% plot data
        surfm(double(latitude), double(longitude), squeeze(corr_spa(:,:,ff,tt)))

        %%% add mooring
        plotm(gps_site(1),gps_site(2),'xk','markersize',16,'linewidth',3)
    %     title(['Frequency: ' num2str(f1(ff)) ' - ' num2str(f2(ff)) ' Hz'])
        colorbar

        caxis([0 0.5])
        
        print(gcf,['./spatial_corr_result/' sitename '/spatial_corr_' ...
            '-freq_' num2str(f1(ff)) '-' num2str(f2(ff)) ...
            '_' datestr(t_beg_num, 'yyyymmdd') '-' datestr(t_end_num, 'yyyymmdd')]  ...
            ,'-dpng')
        
%         saveas(gcf,['./spatial_corr_result/' sitename '/spatial_corr_' ...
%             '-freq_' num2str(f1(ff)) '-' num2str(f2(ff)) ...
%             '_' datestr(t_beg_num, 'yyyymmdd') '-' datestr(t_end_num, 'yyyymmdd')],'fig')
        
    end
    date_loop(1,tt)=t_beg_num;
    date_loop(2,tt)=t_end_num;
    tt=tt+1;  
end

save spatial_cor_results_new_interp.mat corr_spa date_loop f1 f2 latitude longitude

disp('Done')