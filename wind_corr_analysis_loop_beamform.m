clear all
close all
clc

load ANL_SHRU5.mat
t=timestamp_num_spectro;
SPL_ANL_nobeamform=SPL_ANL;

load ANL_SHRU5_beamform.mat

removeOutliers=0;
high_freq=0;
%% Geography limit
sitename = 'SHRU5';
gps_site = [72+54.4580/60 , -(157+29.2442/60)];  


latlimit=[65 85];
dlon=40;
lonlimit=[gps_site(2)-dlon gps_site(2)+dlon];
lonlimit_ok=[lonlimit(2) lonlimit(1)+360];
%% load env_data
wind=T_ecmwf.W10;
sst=T_ecmwf.sst;
tp=T_ecmwf.tp;
tw=timestamp_num_ecmwf;
Nw=length(tw);

% auxData_ok=[wind, sst, tp];
auxData_ok=wind;
%% Limit to a single frequency band 250-350 Hz and 500-1000
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


SPL_ANL_ave=NaN(size(tw)); 
SPL_ANL_ave_up=NaN(size(tw)); 
SPL_ANL_ave_hor=NaN(size(tw)); 

SPL_ANL_ave2=NaN(size(tw)); 
SPL_ANL_ave_up2=NaN(size(tw)); 
SPL_ANL_ave_hor2=NaN(size(tw)); 



for tt=1:Nw-1
    [~, b1]=min(abs(t-tw(tt)));
    [~, b2]=min(abs(t-tw(tt+1)));
    db=b2-b1;
    if db>10
        i_min=b1-round(db);
        if i_min<1
            i_min=1;
        end
        i_max=min([length(t),b1+round(db)]);
        SPL_ANL_ave(tt)=mean(SPL_ANL_ok(i_min:i_max));
        SPL_ANL_ave_up(tt)=mean(SPL_ANL_up_ok(i_min:i_max));
        SPL_ANL_ave_hor(tt)=mean(SPL_ANL_hor_ok(i_min:i_max));
        SPL_ANL_ave2(tt)=mean(SPL_ANL_ok2(i_min:i_max));
        SPL_ANL_ave_up2(tt)=mean(SPL_ANL_up_ok2(i_min:i_max));
        SPL_ANL_ave_hor2(tt)=mean(SPL_ANL_hor_ok2(i_min:i_max));
    end
end
% 
% figure,
% plot(t,SPL_ANL_ok)
% hold on
% plot(tw, SPL_ANL_ave)
% hold on
% plot(tw, SPL_ANL_ave_up)
% hold on
% plot(tw, SPL_ANL_ave_hor)
% datetick
% legend('Original data', 'Averaged data', 'Beamforming up', 'Beamforming horizontally')
% grid on
% title('250 - 350 Hz')
% 
% figure,
% plot(t,SPL_ANL_ok2)
% hold on
% plot(tw, SPL_ANL_ave2)
% hold on
% plot(tw, SPL_ANL_ave_up2)
% hold on
% plot(tw, SPL_ANL_ave_hor2)
% datetick
% legend('Original data', 'Averaged data', 'Beamforming up', 'Beamforming horizontally')
% grid on     
% title('500 - 1000 Hz')    
    
%% Must loop from Nov 1st to May 31

% t_beg=[2016,12,01,00,00,00];  %%% date vector [YY, MM, DD, h, min, sec]
% t_end=[2017,02,01,00,00,00];   %%% date vector [YY, MM, DD, h, min, sec]

t0=[2016,11,01,00,00,00];
t0_num=datenum(t0);

t1=[2017,05,31,00,00,00];
t1_num=datenum(t1);

Nvar_ok=size(auxData_ok,2);
Nloop=length(t0_num:15:t1_num);
date_loop=zeros(2,Nloop);


tt=1;
for t_num_loop=t0_num:15:t1_num
    t_beg_num=t_num_loop;
    t_end_num=t_num_loop+60;

    disp(datestr(t_beg_num, 'yyyymmdd'))



    ind_t_ok=find(tw > t_beg_num & tw < t_end_num);

    auxData_tt=auxData_ok(ind_t_ok,:);      
    vecSPL_tt=SPL_ANL_ave(ind_t_ok);
    vecSPL_tt_up=SPL_ANL_ave_up(ind_t_ok);
    vecSPL_tt_hor=SPL_ANL_ave_hor(ind_t_ok);
    vecSPL_tt2=SPL_ANL_ave2(ind_t_ok);
    vecSPL_tt_up2=SPL_ANL_ave_up2(ind_t_ok);
    vecSPL_tt_hor2=SPL_ANL_ave_hor2(ind_t_ok);
    %%%% Compute correlation

    for ii=1:Nvar_ok
     
        ind_no_nan=(~isnan(auxData_tt(:,ii)) | ~isnan(vecSPL_tt)) ;
        
        %%%% remove corr data which are nearly constant and/or nearly all
        %%%% NaN
        if length(unique(auxData_tt(ind_no_nan,ii)))>10
          
            figure('visible','off');
            h = get(0,'children');
            scrsz = get(0,'ScreenSize');
            set(h,'Position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)])  
            
    
            
%% First frequency band            
            
            ind_plot=1;
            X = [vecSPL_tt(ind_no_nan), auxData_tt(ind_no_nan,ii)];
            X_norma=zscore(X,[],1);            
            tbl = table(X_norma(:,1),X_norma(:,2),'VariableNames', {'X1','X2'});            
            if removeOutliers
                mdl = fitlm(tbl,'X2 ~ X1');
                outliers = mdl.Diagnostics.CooksDistance > 4*mean(mdl.Diagnostics.CooksDistance);
                mdl2 = fitlm(tbl,'X2 ~ X1','Exclude', outliers );       
                [R,pval] = corr(X(~outliers,1),X(~outliers,2),'type','Pearson','rows','all','tail','both');
            else
                mdl2 = fitlm(tbl,'X2 ~ X1');
                [R,pval] = corr(X(:,1),X(:,2),'type','Pearson','rows','all','tail','both');
            end   
                      
            b=subplot(2,3,ind_plot);
            plot(mdl2, 'linewidth', 1);
            xlabel('SPL [250 350] Hz')
            ylabel('Wind speed')
            title('No beamforming')
            legend(b,'off');
            
            plotPos = get(b,'Position');
            if pval < 0.05
                annotation('textbox',plotPos,...
                   'String',strcat(num2str(R,'%3.2f')),...
                   'FontWeight','Bold',...
                   'EdgeColor','none','Tag','corrCoefs','fontsize',14)
            else
                annotation('textbox',plotPos,...
                       'String',strcat(num2str(R,'%3.2f')),...
                       'EdgeColor','none','Tag','corrCoefs','fontsize',14) 
            end
            

            
            ind_plot=2;
            X = [vecSPL_tt_hor(ind_no_nan), auxData_tt(ind_no_nan,ii)];
            X_norma=zscore(X,[],1);            
            tbl = table(X_norma(:,1),X_norma(:,2),'VariableNames', {'X1','X2'});            
            if removeOutliers
                mdl = fitlm(tbl,'X2 ~ X1');
                outliers = mdl.Diagnostics.CooksDistance > 4*mean(mdl.Diagnostics.CooksDistance);
                mdl2 = fitlm(tbl,'X2 ~ X1','Exclude', outliers );       
                [R,pval] = corr(X(~outliers,1),X(~outliers,2),'type','Pearson','rows','all','tail','both');
            else
                mdl2 = fitlm(tbl,'X2 ~ X1');
                [R,pval] = corr(X(:,1),X(:,2),'type','Pearson','rows','all','tail','both');
            end   
                      
            b=subplot(2,3,ind_plot);
            plot(mdl2, 'linewidth', 1);
            xlabel('SPL [250 350] Hz')
            ylabel('Wind speed')
            title('Horizontal beamforming')
            legend(b,'off');
            
            plotPos = get(b,'Position');
            if pval < 0.05
                annotation('textbox',plotPos,...
                   'String',strcat(num2str(R,'%3.2f')),...
                   'FontWeight','Bold',...
                   'EdgeColor','none','Tag','corrCoefs','fontsize',14)
            else
                annotation('textbox',plotPos,...
                       'String',strcat(num2str(R,'%3.2f')),...
                       'EdgeColor','none','Tag','corrCoefs','fontsize',14) 
            end            
            
            
            
            
            
            ind_plot=3;
            X = [vecSPL_tt_up(ind_no_nan), auxData_tt(ind_no_nan,ii)];
            X_norma=zscore(X,[],1);            
            tbl = table(X_norma(:,1),X_norma(:,2),'VariableNames', {'X1','X2'});            
            if removeOutliers
                mdl = fitlm(tbl,'X2 ~ X1');
                outliers = mdl.Diagnostics.CooksDistance > 4*mean(mdl.Diagnostics.CooksDistance);
                mdl2 = fitlm(tbl,'X2 ~ X1','Exclude', outliers );       
                [R,pval] = corr(X(~outliers,1),X(~outliers,2),'type','Pearson','rows','all','tail','both');
            else
                mdl2 = fitlm(tbl,'X2 ~ X1');
                [R,pval] = corr(X(:,1),X(:,2),'type','Pearson','rows','all','tail','both');
            end   
                      
            b=subplot(2,3,ind_plot);
            plot(mdl2, 'linewidth', 1);
            xlabel('SPL [250 350] Hz')
            ylabel('Wind speed')
            title('Vertical (up) beamforming')
            legend(b,'off');
            
            plotPos = get(b,'Position');
            if pval < 0.05
                annotation('textbox',plotPos,...
                   'String',strcat(num2str(R,'%3.2f')),...
                   'FontWeight','Bold',...
                   'EdgeColor','none','Tag','corrCoefs','fontsize',14)
            else
                annotation('textbox',plotPos,...
                       'String',strcat(num2str(R,'%3.2f')),...
                       'EdgeColor','none','Tag','corrCoefs','fontsize',14) 
            end            

%% time series
            if ~high_freq
                ind_plot=4;
                b=subplot(2,3,ind_plot);
                mu=nanmean(vecSPL_tt);
                sigma=nanstd(vecSPL_tt);
                mu_w=nanmean(auxData_tt);
                sigma_w=nanstd(auxData_tt);
                plot(tw(ind_t_ok)-tw(ind_t_ok(1)), (vecSPL_tt-mu)/sigma)
                hold on
                plot(tw(ind_t_ok)-tw(ind_t_ok(1)), (auxData_tt-mu_w)/sigma_w)   
                grid on
                legend('SPL', 'Wind speed')
                grid on
                xlabel('Time (days)')
                
                
                ind_plot=5;
                b=subplot(2,3,ind_plot);
                mu=nanmean(vecSPL_tt_hor);
                sigma=nanstd(vecSPL_tt_hor);
                plot(tw(ind_t_ok)-tw(ind_t_ok(1)), (vecSPL_tt_hor-mu)/sigma)
                hold on
                plot(tw(ind_t_ok)-tw(ind_t_ok(1)), (auxData_tt-mu_w)/sigma_w)   
                grid on
                legend('SPL', 'Wind speed')
                grid on
                xlabel('Time (days)')
                
                
                ind_plot=6;
                b=subplot(2,3,ind_plot);
                mu=nanmean(vecSPL_tt_up);
                sigma=nanstd(vecSPL_tt_up);
                plot(tw(ind_t_ok)-tw(ind_t_ok(1)), (vecSPL_tt_up-mu)/sigma)
                hold on
                plot(tw(ind_t_ok)-tw(ind_t_ok(1)), (auxData_tt-mu_w)/sigma_w)   
                grid on
                legend('SPL', 'Wind speed')
                grid on                
                xlabel('Time (days)')
                
                
            end
            
            
%% second frequency band            
            if high_freq            
            ind_plot=4;
            X = [vecSPL_tt2(ind_no_nan), auxData_tt(ind_no_nan,ii)];
            X_norma=zscore(X,[],1);            
            tbl = table(X_norma(:,1),X_norma(:,2),'VariableNames', {'X1','X2'});            
            if removeOutliers
                mdl = fitlm(tbl,'X2 ~ X1');
                outliers = mdl.Diagnostics.CooksDistance > 4*mean(mdl.Diagnostics.CooksDistance);
                mdl2 = fitlm(tbl,'X2 ~ X1','Exclude', outliers );       
                [R,pval] = corr(X(~outliers,1),X(~outliers,2),'type','Pearson','rows','all','tail','both');
            else
                mdl2 = fitlm(tbl,'X2 ~ X1');
                [R,pval] = corr(X(:,1),X(:,2),'type','Pearson','rows','all','tail','both');
            end   
                      
            b=subplot(2,3,ind_plot);
            plot(mdl2, 'linewidth', 1);
            xlabel('SPL [500 1000] Hz')
            ylabel('Wind speed')
            title('No beamforming')
            legend(b,'off');
            
            plotPos = get(b,'Position');
            if pval < 0.05
                annotation('textbox',plotPos,...
                   'String',strcat(num2str(R,'%3.2f')),...
                   'FontWeight','Bold',...
                   'EdgeColor','none','Tag','corrCoefs','fontsize',14)
            else
                annotation('textbox',plotPos,...
                       'String',strcat(num2str(R,'%3.2f')),...
                       'EdgeColor','none','Tag','corrCoefs','fontsize',14) 
            end
            

            
            ind_plot=5;
            X = [vecSPL_tt_hor2(ind_no_nan), auxData_tt(ind_no_nan,ii)];
            X_norma=zscore(X,[],1);            
            tbl = table(X_norma(:,1),X_norma(:,2),'VariableNames', {'X1','X2'});            
            if removeOutliers
                mdl = fitlm(tbl,'X2 ~ X1');
                outliers = mdl.Diagnostics.CooksDistance > 4*mean(mdl.Diagnostics.CooksDistance);
                mdl2 = fitlm(tbl,'X2 ~ X1','Exclude', outliers );       
                [R,pval] = corr(X(~outliers,1),X(~outliers,2),'type','Pearson','rows','all','tail','both');
            else
                mdl2 = fitlm(tbl,'X2 ~ X1');
                [R,pval] = corr(X(:,1),X(:,2),'type','Pearson','rows','all','tail','both');
            end   
                      
            b=subplot(2,3,ind_plot);
            plot(mdl2, 'linewidth', 1);
            xlabel('SPL [500 1000] Hz')
            ylabel('Wind speed')
            title('Horizontal beamforming')
            legend(b,'off');
            
            plotPos = get(b,'Position');
            if pval < 0.05
                annotation('textbox',plotPos,...
                   'String',strcat(num2str(R,'%3.2f')),...
                   'FontWeight','Bold',...
                   'EdgeColor','none','Tag','corrCoefs','fontsize',14)
            else
                annotation('textbox',plotPos,...
                       'String',strcat(num2str(R,'%3.2f')),...
                       'EdgeColor','none','Tag','corrCoefs','fontsize',14) 
            end            
                        
            
            
            
            ind_plot=6;
            X = [vecSPL_tt_up2(ind_no_nan), auxData_tt(ind_no_nan,ii)];
            X_norma=zscore(X,[],1);            
            tbl = table(X_norma(:,1),X_norma(:,2),'VariableNames', {'X1','X2'});            
            if removeOutliers
                mdl = fitlm(tbl,'X2 ~ X1');
                outliers = mdl.Diagnostics.CooksDistance > 4*mean(mdl.Diagnostics.CooksDistance);
                mdl2 = fitlm(tbl,'X2 ~ X1','Exclude', outliers );       
                [R,pval] = corr(X(~outliers,1),X(~outliers,2),'type','Pearson','rows','all','tail','both');
            else
                mdl2 = fitlm(tbl,'X2 ~ X1');
                [R,pval] = corr(X(:,1),X(:,2),'type','Pearson','rows','all','tail','both');
            end   
                      
            b=subplot(2,3,ind_plot);
            plot(mdl2, 'linewidth', 1);
            xlabel('SPL [500 1000] Hz')
            ylabel('Wind speed')
            title('Vertical (up) beamforming')
            legend(b,'off');
            
            plotPos = get(b,'Position');
            if pval < 0.05
                annotation('textbox',plotPos,...
                   'String',strcat(num2str(R,'%3.2f')),...
                   'FontWeight','Bold',...
                   'EdgeColor','none','Tag','corrCoefs','fontsize',14)
            else
                annotation('textbox',plotPos,...
                       'String',strcat(num2str(R,'%3.2f')),...
                       'EdgeColor','none','Tag','corrCoefs','fontsize',14) 
            end            
            end
            
            

        end
        
                               
    end    
    
    

    print(gcf,['./local_corr_result/' sitename '_beamform/wind_corr' ...
        '_' datestr(t_beg_num, 'yyyymmdd') '-' datestr(t_end_num, 'yyyymmdd')]  ...
        ,'-dpng') 
    
    date_loop(1,tt)=t_beg_num;
    date_loop(2,tt)=t_end_num;
    tt=tt+1;  
end

% save spatial_cor_results_beamform.mat corr_spa_ave corr_spa_ave_up corr_spa_ave_hor date_loop f1 f2 latitude longitude SPL_ANL_ave SPL_ANL_ave_up SPL_ANL_ave_hor

disp('Done')