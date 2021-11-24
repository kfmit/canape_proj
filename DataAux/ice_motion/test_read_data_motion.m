close all
clear all
clc
fileID = fopen('./data_nsidc/icemotion.grid.daily.2016001.n.v3.bin','r','ieee-le');
% fileID = fopen('icemotion.grid.month.2017.01.n.v3.bin','r','ieee-le');
toto = fread(fileID, Inf,'int16');
fclose(fileID);


u=toto(1:3:end)/10;
v=toto(2:3:end)/10;
err=toto(3:3:end);

%%% Arctic data are 361*361
N=361;
u_ok=reshape(u,[N,N]);
v_ok=reshape(v,[N,N]);


[x,y]=meshgrid(1:N,1:N);

figure
quiver(x,y,u_ok,v_ok)
xlim([1 N])
ylim([1 N])
grid on

%% Coordinate
fileID = fopen('north_x_y_lat_lon.txt');
coord = fscanf(fileID, '%d %d %f %f',[4, Inf]);
fclose(fileID);


lat=reshape(coord(3,:,:),[N,N]);
lon=reshape(coord(4,:,:),[N,N]);


% lat=fliplr(lat);
% lon=fliplr(lon);

%%%%% Same as reshape
% coord(1:2,:,:)=coord(1:2,:,:)+1; %%% In matlab, index starts at 1
% lat=zeros(N,N);
% lon=zeros(N,N);
% for ii=1:size(coord,2)
%     ligne=coord(1,ii);
%     col=coord(2,ii);
%     lat(ligne,col)=coord(3,ii);
%     lon(ligne,col)=coord(4,ii);    
% end

%% Map
fileID = fopen('nsidc_north_map.txt');

map0 = fread(fileID, Inf);
fclose(fileID);

%%% Arctic map is 722*722
Nmap=722;
map=reshape(map0, [Nmap,Nmap]);
map=fliplr(map); %%% must turn map to fit it to data (cf show_vectors_v3.pro)

map(map==184)=0;
map(map==255)=1;

map_scale=linspace(0,1,Nmap)*(N-1)+1;

figure
imagesc(map_scale, map_scale,map)
colormap gray


%% Map and data
figure
imagesc(map_scale, map_scale,map)
colormap gray
hold on
quiver(x,y,u_ok,v_ok)
axis xy

% figure
% subplot(131)
% imagesc(lat)
% title('latitude')
% axis xy
% colorbar
% subplot(132)
% imagesc(lon)
% title('longitude')
% axis xy
% colorbar
% subplot(133)
% imagesc(map_scale, map_scale,map)
% hold on
% quiver(x,y,u_ok,v_ok)
% axis xy
