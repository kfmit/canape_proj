%% Geographic Location of SHRUs

clear
clc
close all

% LATITUDES
shrulat = [54.4123 45.2347 40.6924 36.6582 45.4580]/60;         % N
shrulat = shrulat+72;
shrulon =[159.018066666666667 158.2721 157.9108 157.5375 157.4874 ];    % W
shrudepth = [302 308 303 301 445]*-1;
shrutime = [21:04:00 23:37:00 02:04:00 04:21:00 23:28:00];
shrudate = ['24-Oct-2016' '24-Oct-2016' '25-Oct-2016' '25-Oct-2016' '21-Oct-2016'];


SHRUinfo = struct('lat',shrulat,'lon',shrulon,'depth',shrudepth,'time',shrutime,'date',shrudate);

% PLOT
figure(1)
geoscatter(SHRUinfo.lat(1),SHRUinfo.lon(1),100,'^','filled')
hold on
geoscatter(SHRUinfo.lat(2),SHRUinfo.lon(2),100,'^','filled')
hold on
geoscatter(SHRUinfo.lat(3),SHRUinfo.lon(3),100,'^','filled')
geoscatter(SHRUinfo.lat(4),SHRUinfo.lon(4),100,'^','filled')
geoscatter(SHRUinfo.lat(5),SHRUinfo.lon(5),100,'^','filled')
hold off
legend('SHRU1','SHRU2', 'SHRU3', 'SHRU4', 'SHRU5')
title('Location of Hydrophone Arrays')
gb.SizeLegendTitle = 'Maximum Height';
geobasemap colorterrain

%% Resolution of Whole Array


figure(2)
scatter3(SHRUinfo.lat(1),SHRUinfo.lon(1),SHRUinfo.depth(1),100,'^','filled')
hold on
scatter3(SHRUinfo.lat(2),SHRUinfo.lon(2),SHRUinfo.depth(2),100,'^','filled')
scatter3(SHRUinfo.lat(3),SHRUinfo.lon(3),SHRUinfo.depth(3),100,'^','filled')
scatter3(SHRUinfo.lat(4),SHRUinfo.lon(4),SHRUinfo.depth(4),100,'^','filled')
scatter3(SHRUinfo.lat(5),SHRUinfo.lon(5),SHRUinfo.depth(5),100,'^','filled')
hold off

zlim([-500 0])
xlabel('Latitude ^\circ N')
ylabel('Longitude ^\circ W')
zlabel('Depth (m)')
title('Hydrophone Location and Depth')
