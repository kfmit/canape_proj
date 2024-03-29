%% MUST RUN ICE COVERAGE CORR FIRST
% run kat_ice_coverage_corr_cutter.m FIRST

close all
figure(2)
t=tiledlayout(3,3)

freq_array1=[275 475 975 1475];
freq_array2=[325 525 1025 1525];

loop_end = 11;

for tt=3:loop_end
    %% CUTTING VALUES
    %%% plot data over 0.5
    cut_val = 0.40
%     corr_spa_ave2_cut=zeros(119,177,tt);

    t_beg_num=date_loop(1,tt);
    t_end_num=date_loop(2,tt);
    % take out anything below cutval
    for i = 1:119
        for ii = 1:177
            if 1==isnan(corr_spa_ave2(i,ii,tt))
                corr_spa_ave2_cut(i,ii,tt)=nan;
                % do nothing
            elseif (corr_spa_ave2(i,ii,tt))<=cut_val
                % change tht index to a nan if less than 0.4
                corr_spa_ave2_cut(i,ii,tt)=nan;
            elseif (corr_spa_ave2(i,ii,tt))>cut_val
                corr_spa_ave2_cut(i,ii,tt)=corr_spa_ave2(i,ii,tt);
            end
        end
    end

    % clear so not rewritten every round
end
%%

    %%%%%%%%%%%%%% INSERT ISLAND REMOVER HERE %%%%%%%%%%%%%%
corr_spa_ave2_cut(40,70,3) = nan;
corr_spa_ave2_cut(27,90,7) = nan;
corr_spa_ave2_cut(27,90,8) = nan;
corr_spa_ave2_cut(32,80,9) = nan;
corr_spa_ave2_cut(27,90,11) = nan;
corr_spa_ave2_cut(29,85,11) = nan;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
    for tt = 3:11
    clear r_ind c_ind
    clear dist_map
    clear nans_corr_cut 


        nans_corr_cut = ~isnan(corr_spa_ave2_cut(:,:,tt)); % nan is 0, # is 1
        [r_ind, c_ind,v_corr] = find(nans_corr_cut);    % finds non zeros (ie the #s) and saves inds
        % length x whatever lat/long area is shows coverage
        uncut_corr = corr_spa_ave2_cut(r_ind,c_ind,tt);   %

        % find size of boi
        earthellipsoid = referenceSphere('earth','km');
        %     area_base = areaquad(double(latitude(1,1)),double(longitude(1,1)),double(latitude(1,2)),double(longitude(1,2)),earthellipsoid,'degrees');
        %     area_cov(tt) = areaint(double(latitude(r_ind,c_ind)),double(longitude(r_ind,c_ind)),earthellipsoid)

        % count # of pixels active, mulitply by pixel
        %     r_ind*area
        toto=find(corr_spa_ave2_cut(:,:,tt)>0.4) ;
        N_pixel=length(toto);
        pixel_lat=latitude(toto);
        pixel_lon=longitude(toto);
        area_cov(tt) = N_pixel*62.5*62.5    % this is in km!!!!


        % find all the distances
        for iv=1:length(r_ind)
            ri=r_ind(iv);
            ci=c_ind(iv);
            dist_map(iv)=distance(gps_site(1), gps_site(2), double(latitude(ri, ci)),double(longitude(ri,ci)),referenceSphere('Earth'));
        end

        % find index of closest point
        [min_dist(tt), min_ind(tt)]=min(dist_map);
        minlat(tt)=double(latitude(r_ind(min_ind(tt)), c_ind(min_ind(tt))));
        minlon(tt)=double(longitude(r_ind(min_ind(tt)),c_ind(min_ind(tt))));

        % find index of farthest
        [max_dist(tt),max_ind(tt)]=max(dist_map);
        maxlat(tt)=double(latitude(r_ind(max_ind(tt)), c_ind(max_ind(tt))));
        maxlon(tt)=double(longitude(r_ind(max_ind(tt)),c_ind(max_ind(tt))));

        % find index of mid of avg point ( for plotting purposes)
        avg_dist(tt)=mean(dist_map);
        [close_avg_dist(tt), avg_ind(tt)]=min(abs(mean(dist_map)-dist_map));
        avglat(tt)=double(latitude(r_ind(avg_ind(tt)), c_ind(avg_ind(tt))));
        avglon(tt)=double(longitude(r_ind(avg_ind(tt)),c_ind(avg_ind(tt))));

        if freq_range1 == 1475
            avglon(7)=avglon2(7);
            avglat(7)=avglat2(7);
        end
    end

    %% MEGAMAP MAKER %%%%%%%%%%%%%%%%%%%%%%%%%5
    %%%% make the figure(s) plotting this
    for tt = 3:11
        t_beg_num=date_loop(1,tt);
        t_end_num=date_loop(2,tt);
        %     figure('visible','off');
        fig_name= figure(2)
        nexttile
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

        %%% plot CUTOFF data
        surfm(double(latitude), double(longitude), squeeze(corr_spa_ave2_cut(:,:,tt)))

        %%% add mooring
        plotm(gps_site(1),gps_site(2),'xk','markersize',16,'linewidth',3)
        [toto, tata]=ind2sub(size(latitude), indc_ave2(tt));
        %%% add location of max corr point
        plotm(double(latitude(toto, tata)),double(longitude(toto,tata)),'xr','markersize',10,'linewidth',3)

        %%% location of far point
        px = plotm(maxlat(tt),maxlon(tt),'.m','markersize',20,'linewidth',3)

        %%% location of min point
        pn = plotm(minlat(tt),minlon(tt),'.c','markersize',20,'linewidth',3)

        %%% location of avg point
        %     pa = plotm(avglat(tt),avglon(tt),'dg','markersize',10,'linewidth',3)%'MarkerFaceColor','#D95319')
        pa = plotm(avglat(tt),avglon(tt),'dg','markersize',10,'linewidth',3)%'MarkerFaceColor','#D95319')

        %     %%% add edges
        %     plotm(latitude_edge_ok, longitude_edge_ok, '.k','markersize',3,'linewidth',1)


        title([datestr(t_beg_num, 'dd mmm yyyy') ' to ' datestr(t_end_num, 'dd mmm yyyy')])
        xlabel(['Min: ' num2str(min_dist(tt)) ' Max: ' num2str(max_dist(tt))]) %' Avg: ' num2str(avg_dist(tt))])

        sgtitle(['Correlation over ' num2str(cut_val) ' for ' num2str(freq_range1) '-' num2str(freq_range2) ' Hz'])

    end

    c=colorbar;
    c.Label.String = 'Correlation coefficient';
    caxis([0.4 0.7]);
    c.Layout.Tile = 'east';
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
    ll = legend([pn px pa],{'Minimum', 'Maximum', 'Average'})
    ll.Layout.Tile = 'north'
    ll.Title.String = 'Distance Markers'
    ll.FontSize = 10;

    %     scrsz = get(0,'ScreenSize');
    % set(fig_name,'Position',scrsz)
    %         print(gcf,['./new_figs/spatial_corr_result_50/megamap_' num2str(freq_range1) '_' ...
    %             num2str(freq_range2) '_' num2str(cut_val)],'-dpng')

    %%%% end of the figure generator %%%%
    % savestring = ['cutoff_' num2str(freq_range1) '_' num2str(freq_range2) '_icecorr']
    %     save(savestring,'freq_range1','freq_range2','corr_spa_ave2_cut','area_cov','min_dist','minlat','minlon',...
    %         'max_dist','maxlat','maxlon','avg_dist','avglat','avglon')


    %%
    %%%%% Functions! %%%%%%%
    dist_shrus=distance(gps_site(1), gps_site(2),gps_site_shru1(1), gps_site_shru1(2),referenceSphere('Earth'))/1000;

    dist_corr_shru5=[min(dist(3:loop_end)) mean(dist(3:loop_end)) max(dist(3:loop_end))]/1000;

    % dist_corr_shru1=[min(dist_shru1(3:loop_end)) mean(dist_shru1(3:loop_end)) max(dist_shru1(3:loop_end))]/1000

    r_shru5=[min(R(3:loop_end)) mean(R(3:loop_end)) max(R(3:loop_end))];

    % r_shru1=[min(R_shru1(3:loop_end)) mean(R_shru1(3:loop_end)) max(R_shru1(3:loop_end))]

