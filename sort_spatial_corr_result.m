close all
clear all
clc


shru_sn='SHRU5';
fig_name='SpatialPearson_chn_0_freq_50-500_p15';
dest_folder='./spatial_corr_result';

toto=dir(['./Results/**/Figures/' fig_name '.png']);
N=length(toto);

ind_ok=[];
for nn=1:N
    if contains(toto(nn).folder,shru_sn)
        ind_ok=[ind_ok nn];
    end
end

N_ok=length(ind_ok);

% dest_folder

for nn=1:N_ok
    file=toto(ind_ok(nn)).name;
    folder=toto(ind_ok(nn)).folder;
    
    copyfile([folder '/' file], [dest_folder '/'  shru_sn '/' fig_name '_' folder(end-24:end-8) '.png'] );
end
    
