clc
clear all
close all

NameDataset = strcat('WHOI_Bonnel_201706.nc');

info=ncinfo(NameDataset);

varNames = {info.Variables.Name}; % ne fonctionne pas

varNames = {'time','longitude','latitude',...
        'mwp','swh','u10','v10','sst','tp'};% select only a few variables

for nnn=1:length(varNames)
    eval([varNames{nnn} '=ncread(NameDataset,''' varNames{nnn} ''');'])
end

time_str = double(time)/24 + datenum('1900-01-01 00:00:00');
time_str = sort(time_str,'ascend');

ind = find(~cellfun(@isempty,(strfind(varNames,'time'))),1,'first'); 
varNames(ind) = [];
ind = find(~cellfun(@isempty,(strfind(varNames,'longitude'))),1,'first'); 
varNames(ind) = [];
ind = find(~cellfun(@isempty,(strfind(varNames,'latitude'))),1,'first'); 
varNames(ind) = [];

datestr(time_str(:),'yyyy-mm-dd-HH') 
disp(['Your ECMWF data are from ' datestr(time_str(1),'yyyy-mm-dd-HH') ' to ' datestr(time_str(end),'yyyy-mm-dd-HH') ])
pause

cpt=1;
for tt=1:length(time_str)
tt/length(time_str)
    
    NameFile = datestr(time_str(tt,:),'yyyy-mm-dd-HH');
    
    for nnn=1:length(varNames)
        eval(['str.' varNames{nnn} ' = ' varNames{nnn} '(:,:,tt);'])
    end
    
    str.latitude = latitude;
    str.longitude = longitude;

%     save(strcat(cd,filesep,NameFile,'.mat'),'-struct','str','-v7.3')
end
