close all
clear all
clc

load spatial_cor_results_interp_new.mat
gps_site = [72+54.4580/60 , -(157+29.2442/60)];  
dlon=40;
lonlimit=[gps_site(2)-dlon gps_site(2)+dlon];
lonlimit_ok=[lonlimit(2) lonlimit(1)+360];
latlimit=[65 85];

load ANL_SHRU5.mat
t=timestamp_num_spectro;

% load auxData_icedrift_on_psd_time.mat

load var_osisaf.mat
t_osisaf=datenum_osisaf;

%% Reject area of no-interest
toto=find(latitude<=latlimit(1) | latitude>=latlimit(2));
tata=find(longitude >= lonlimit_ok(1) & longitude <=lonlimit_ok(2));
tutu=union(toto, tata);
%%%%% verification
% lat_ok=latitude;
% lon_ok=longitude;
% lat_ok(tutu)=NaN;
% lon_ok(tutu)=NaN;
% figure, imagesc(lat_ok)
% figure, imagesc(lon_ok)



% %% Interpolated Ice data
Nloop=size(date_loop,2);
ff=2;
% 
% figure
% for tt=1:Nloop
%    J=squeeze(corr_spa_near(:,:,tt));
%    J(tutu)=NaN;
% %    figure, imagesc(J)  
%    [cmax(tt), ind]=max(abs(J(:)));
%    [ii jj]=ind2sub(size(latitude), ind);
%    
%    
%    t_beg_num=date_loop(1,tt);
%    t_end_num=date_loop(2,tt);   
%    ind_t_ok=find(t >= t_beg_num & t <= t_end_num);     
%    ind_t_ok_osisaf=find(t_osisaf >= t_beg_num & t_osisaf <= t_end_num);
%   
%    SPL_ANL_ok=SPL_ANL(ind_t_ok,ff);
%    d_ok=squeeze(d(ind_t_ok_osisaf,ii,jj));   
% %    d_ok_interp=interp1(t_osisaf(ind_t_ok_osisaf), d_ok, t(ind_t_ok), 'nearest');
%    d_ok_interp=interp1(t_osisaf(ind_t_ok_osisaf), d_ok, t(ind_t_ok));
%    t_ok=t(ind_t_ok); %%% common time axis
%    
%       
%    ind_no_nan=~isnan(d_ok_interp);
%    
% 
%     [SPL_ANL_ok_norm, mu_spl, sigma_spl]=zscore(SPL_ANL_ok(ind_no_nan));
%     [d_ok_norm, mu_d, sigma_d]=zscore(d_ok_interp(ind_no_nan));
%     
%     [R(tt), Pvalue]=corr(SPL_ANL_ok_norm,d_ok_norm,'type','Pearson','rows','all','tail','both');
%     
% %     figure
%     subplot(5,3,tt)
%     yyaxis left
%     plot(t_ok,(SPL_ANL_ok-mu_spl)/sigma_spl)
%     yyaxis right
%     plot(t_ok(ind_no_nan), d_ok_norm, '.')
%     xlim([t_ok(1) t_ok(end)])
% %     datetick2('x', 'mm/dd')
%     grid on
%     title(['Correlation ' num2str(cmax(tt)) ' / ' num2str(R(tt))])
%    
% end

%% Averaged SPL




figure
for tt=1:Nloop
   J=squeeze(corr_spa_ave2(:,:,tt));
   J(tutu)=NaN;
%    figure, imagesc(J)  
   [cmax(tt), ind]=max(abs(J(:)));
   [ii jj]=ind2sub(size(latitude), ind);
   
   
   t_beg_num=date_loop(1,tt);
   t_end_num=date_loop(2,tt);      
   ind_t_ok_osisaf=find(t_osisaf >= t_beg_num & t_osisaf <= t_end_num);
  
   SPL_ANL_ok=SPL_ANL_ave2(ind_t_ok_osisaf);
   d_ok=squeeze(d(ind_t_ok_osisaf,ii,jj));   
%    d_ok_interp=interp1(t_osisaf(ind_t_ok_osisaf), d_ok, t(ind_t_ok), 'nearest');
   t_ok=t_osisaf(ind_t_ok_osisaf); %%% common time axis
   
      
    ind_no_nan=~isnan(d_ok);
   

    [SPL_ANL_ok_norm, mu_spl, sigma_spl]=zscore(SPL_ANL_ok(ind_no_nan));
    [d_ok_norm, mu_d, sigma_d]=zscore(d_ok(ind_no_nan));
    d_ok_norm2=(d_ok-mu_d)/sigma_d;
    
    [R(tt), Pvalue]=corr(SPL_ANL_ok_norm,d_ok_norm,'type','Pearson','rows','all','tail','both');
    
    toto=find(~isnan(d_ok));
    
    Ntoto=length(toto);
    nb=1;
    block{nb}=d_ok_norm2(toto(1));
    ind_b{nb}=toto(1);
    for iii=2:Ntoto
        if toto(iii)-toto(iii-1) == 1
            block{nb}=[block{nb} d_ok_norm2(toto(iii))];
            ind_b{nb}=[ind_b{nb} toto(iii)];
        else %%% new block
            nb=nb+1;
            block{nb}= d_ok_norm2(toto(iii));     
            ind_b{nb}= toto(iii);
        end
    end
                
    
%     figure
    subplot(5,3,tt)
    plot(t_ok-t_ok(1),(SPL_ANL_ok-mu_spl)/sigma_spl, 'Color', [0    0.4470    0.7410])
       
    hold on
    for bb=1:nb
        plot(t_ok(ind_b{bb})-t_ok(1), block{bb}, '-', 'Color', [0.8500    0.3250    0.0980])
        hold on
    end
    
    plot(t_ok(ind_no_nan)-t_ok(1), d_ok_norm, 'o', 'Color', [0.8500    0.3250    0.0980])
    grid on
    title({[datestr(t_ok(1), 'dd mmm yyyy') ' to ' datestr(t_ok(end), 'dd mmm yyyy')] , ['R=' num2str(R(tt))]})
end

