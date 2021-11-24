
clc
clear all
close all

dirnc = dir('./*.nc');

% sitename = 'SHRU1';
% gps_site = [72+54.4123/60 , -(159+1.0840/60)];

% sitename = 'SHRU2';
% gps_site = [72+45.2347/60 , -(158+16.3243/60)];

% sitename = 'SHRU3';
% gps_site = [72+40.6924/60 , -(157+54.6493/60)];
%
% sitename = 'SHRU4';
% gps_site = [72+36.6582/60 , -(157+32.2475/60)];
%
sitename = 'SHRU5';
gps_site = [72+54.4580/60 , -(157+29.2442/60)];


visu_map = 0;



all_T=[];
for dd = 1:size(dirnc,1)

    NameDataset = ['./' dirnc(dd).name];

    info=ncinfo(NameDataset);
    
    varNames = {info.Variables.Name};% take all variables available
       
    varNames = {'time','longitude','latitude',...
        'u10','v10','sst','tp'};% select only a few variables

    for nnn=1:length(varNames)
        nnn/length(varNames)
        eval([varNames{nnn} '=ncread(NameDataset,''' varNames{nnn} ''');'])
    end

    time_str = double(time)/24 + datenum('1900-01-01 00:00:00');
    time_str = sort(time_str,'ascend');

    %%% on vire les colonnes time, longitude et latitude
    ind = find(~cellfun(@isempty,(strfind(varNames,'time'))),1,'first'); 
    varNames(ind) = [];
    ind = find(~cellfun(@isempty,(strfind(varNames,'longitude'))),1,'first'); 
    varNames(ind) = [];
    ind = find(~cellfun(@isempty,(strfind(varNames,'latitude'))),1,'first'); 
    varNames(ind) = [];

    datestr(time_str(:),'yyyy-mm-dd-HH'); 
    disp(['Your ECMWF data are from ' datestr(time_str(1),'yyyy-mm-dd-HH') ' to ' datestr(time_str(end),'yyyy-mm-dd-HH') ])

    longitude(longitude>180)=longitude(longitude>180)-360;
       
    [c_lat,ind_lat] = min(abs(gps_site(1)-latitude));
    [c_lon,ind_lon] = min(abs(gps_site(2)-longitude));

    var = zeros(length(time_str),length(varNames));
    datevect = {};
    for tt=1:length(time_str)
    tt/length(time_str)

        NameFile = datestr(time_str(tt,:),'yyyymmddHHMMSS');


        if visu_map
            figure, hold on,            
            imagesc(u10(:,:,tt))
            plot(ind_lat,ind_lon,'xr','markersize',16,'linewidth',3)
            pause
            
        end
            
        datevect{tt,1} = NameFile;
    end

    for nnn=1:length(varNames)
       eval(['var(:,nnn)=',varNames{nnn},'(ind_lon,ind_lat,:);' ]);
    end
    T=table(datevect,var(:,1),var(:,2),var(:,3),var(:,4)...
        ,'VariableNames',[{'Time'},varNames(1:end)]);
    
    T(end,:)=[];

    all_T = [all_T ; T];
end

% redefine variables in T
final_T = [all_T(:,1) , table(sqrt( all_T.u10 .^2 + all_T.v10.^2)) , all_T(:,[4,5])];
final_T.Properties.VariableNames = {'timestamp','W10','sst','tp'};

[~,ind]=sort(final_T{:,1});
final_T=final_T(ind,:);

writetable(final_T,[cd '/variables_ECMWF_' sitename '.csv']) 

