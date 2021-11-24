close all
clear all
clc

t1=datenum(2016,10,01);
t2=datenum(2017,08,31);

N=t2-t1;

prefix='wget -P smos_sea_ice_thickness http://icdc.cen.uni-hamburg.de/thredds/fileServer/ftpthredds/smos_sea_ice_thickness/v3/SMOS_Icethickness_v3.1_north_';

fileID = fopen('wget_all_files.txt','w');
for nn=0:N
    t=t1+nn;
    datestring=datestr(t,'yyyymmdd');
    str=strcat(prefix,datestring,'.nc');
    fprintf(fileID, '%s', str);
    fprintf(fileID, '\n');
end