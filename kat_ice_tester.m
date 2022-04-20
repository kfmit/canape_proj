% find size of area
A = [1 2 nan 3]
testA = ~isnan(A)
testind = find(testA)
Aproof = A(testind)

%%
nans_corr_cut = ~isnan(corr_spa_ave2_cut(:,:,3)); %nan is 0, # is 1
[r_ind, c_ind,v_corr] = find(nans_corr_cut);    % finds non zeros and saves inds
% length x whatever lat/long area is shows coverage
uncut_corr = corr_spa_ave2_cut(r_ind,c_ind,3);   %

% find all the distances
for iv=1:length(r_ind)
    ri=r_ind(iv);
    ci=c_ind(iv);
dist_map(iv)=distance(gps_site(1), gps_site(2), double(latitude(ri, ci)),double(longitude(ri,ci)),referenceSphere('Earth'));
end
% find index of closest point
[min_dist, min_ind]=min(dist_map);
% find index of farthest
[max_dist,max_ind]=max(dist_map);
% find index of mid of avg point ( for plotting purposes)
[avg_dist, avg_ind]=min(abs(mean(dist_map)-dist_map));