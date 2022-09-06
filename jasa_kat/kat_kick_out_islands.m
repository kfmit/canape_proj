% test cutter
close all

toto=squeeze(corr_spa_ave2(:,:,11));
toto(toto<0.4)=NaN;
figure
imagesc(squeeze(corr_spa_ave2_cut(:,:,11)))
figure
imagesc(toto)
colorbar
% toto(55,83) = nan;
% toto(53,83) = nan;
% toto(49,84) = nan;
% toto(52,84) = nan;
% toto(54,84) = nan;
% toto(56,84) = nan;
% toto(53,85) = nan;
% toto(54,85) = nan;
toto(32,80) = nan;
figure
imagesc(toto)

% first run kat_ice_coverage_corr_cutter.m for select frequency
% then come here and copy paste in cmd window

%% 300 Hz
corr_spa_ave2_cut(29,85,7) = nan;
corr_spa_ave2_cut(55,83,8) = nan;
corr_spa_ave2_cut(53,83,8) = nan;
corr_spa_ave2_cut(49,84,8) = nan;
corr_spa_ave2_cut(52,84,8) = nan;
corr_spa_ave2_cut(54,84,8) = nan;
corr_spa_ave2_cut(53,85,8) = nan;
corr_spa_ave2_cut(54,85,8) = nan;
corr_spa_ave2_cut(56,85,8) = nan;
corr_spa_ave2_cut(47,90,9) = nan;
corr_spa_ave2_cut(27,90,11) = nan;

%% 500 Hz
corr_spa_ave2_cut(28,87,7) = nan;
corr_spa_ave2_cut(29,85,7) = nan;
corr_spa_ave2_cut(55,83,8) = nan;
corr_spa_ave2_cut(52,84,8) = nan;
corr_spa_ave2_cut(54,84,8) = nan;
corr_spa_ave2_cut(55,84,8) = nan;
corr_spa_ave2_cut(56,84,8) = nan;
corr_spa_ave2_cut(53,85,8) = nan;
corr_spa_ave2_cut(54,85,8) = nan;
corr_spa_ave2_cut(56,85,8) = nan;
corr_spa_ave2_cut(56,86,8) = nan;
corr_spa_ave2_cut(55,82,8) = nan;


%% 1000 Hz
corr_spa_ave2_cut(40,70,3) = nan;
corr_spa_ave2_cut(55,83,8) = nan;
corr_spa_ave2_cut(27,90,8) = nan;
corr_spa_ave2_cut(32,80,8) = nan;
corr_spa_ave2_cut(32,82,11) = nan;
corr_spa_ave2_cut(31,86,11) = nan;
corr_spa_ave2_cut(32,82,11) = nan;
corr_spa_ave2_cut(31,86,11) = nan;
corr_spa_ave2_cut(27,90,11) = nan;


%% 1500 Hz
corr_spa_ave2_cut(40,70,3) = nan;
corr_spa_ave2_cut(27,90,7) = nan;
corr_spa_ave2_cut(27,90,8) = nan;
corr_spa_ave2_cut(32,80,9) = nan;
corr_spa_ave2_cut(27,90,11) = nan;
corr_spa_ave2_cut(29,85,11) = nan;

%% test cutter
% 
% for i_lat = 1:119
%     for i_lon = 1:177
%         if toto()
%             
%         end
%     end
% end