clear all
close all
clc

mooring='1';
z=[22     40     55     70     85    100    115];  %%% depths for SBE56
Nz=length(z);

%% SBE37-1
n_capt=1;
name=['SBE37-UDel' mooring '-1.hourly'];
fid=fopen(name);
datacell=textscan(fid, '%f %f %f %f %f %f', 'HeaderLines', 16);
fclose(fid);

time_jd{n_capt}=datacell{1,1}; %%% Julian-Day, 1.5 is noon Jan.-1 2016 --> 0 is Jan 0 2016 at 00h00
pressure=datacell{1,2}; %%% dbar
temp_pot=datacell{1,3}; %%% potential temperature, degree Celsius
sal{n_capt}=datacell{1,4}; %%% psu
temp{n_capt}=datacell{1,5}; %%% potential temperature, degree Celsius
rho_pot_ano=datacell{1,6}; %%% Potential_Density_Anomaly, kg/m^3

rho_pot{n_capt}=rho_pot_ano+1000; %%% Potential_Density, kg/m^3
clear rho_pot_ano

%%%% NB "potential" means as if meansured at the surface (without pressure
%%%% burden)
g=9.806;
depth{n_capt}=((((-1.82e-15*pressure + 2.279e-10).*pressure-2.2512e-5).*pressure + 9.72659).*pressure)/g; %%% from SBE website

t0=datenum([2016,1,0,0,0,0]); %%% t=0 in Julian Days
time_num{n_capt}=t0+time_jd{n_capt};


%% SBE37-2
n_capt=2;
name=['SBE37-UDel' mooring '-2.hourly'];
fid=fopen(name);
datacell=textscan(fid, '%f %f %f %f %f %f', 'HeaderLines', 16);
fclose(fid);

time_jd{n_capt}=datacell{1,1}; %%% Julian-Day, 1.5 is noon Jan.-1 2016 --> 0 is Jan 0 2016 at 00h00
pressure=datacell{1,2}; %%% dbar
temp_pot=datacell{1,3}; %%% potential temperature, degree Celsius
sal{n_capt}=datacell{1,4}; %%% psu
temp{n_capt}=datacell{1,5}; %%% potential temperature, degree Celsius
rho_pot_ano=datacell{1,6}; %%% Potential_Density_Anomaly, kg/m^3

rho_pot{n_capt}=rho_pot_ano+1000; %%% Potential_Density, kg/m^3
clear rho_pot_ano

%%%% NB "potential" means as if meansured at the surface (without pressure
%%%% burden)
g=9.806;
depth{n_capt}=((((-1.82e-15*pressure + 2.279e-10).*pressure-2.2512e-5).*pressure + 9.72659).*pressure)/g; %%% from SBE website

% t0=datenum([2016,1,0,0,0,0]); %%% same as above
time_num{n_capt}=t0+time_jd{n_capt};

%% SBE56
name=['SBE56-UDel' mooring '.hourly'];
fid=fopen(name);
datacell=textscan(fid, '%f %f %f %f %f %f %f %f', 'HeaderLines', 17);
fclose(fid);

for n_capt=3:9
    time_jd{n_capt}=datacell{1,1}; %%% Julian-Day, 1.5 is noon Jan.-1 2016 --> 0 is Jan 0 2016 at 00h00
%     t0=datenum([2016,1,0,0,0,0]); %%% same as above
    time_num{n_capt}=t0+time_jd{n_capt};
    temp{n_capt}=datacell{1,n_capt-1};
    depth{n_capt}=z(n_capt-2);
end

%% Plot

for nn=1:9
    zzz(nn)=mean(depth{nn});
end
[~, b]=sort(zzz);

toto=[];
figure
for nn=1:9
    ind_capt=b(nn);
    p{nn}=subplot(9,1,nn);
    plot(time_num{ind_capt},temp{ind_capt})
    grid on
    title(['Mooring #' mooring ' - Depth = ' num2str(mean(depth{ind_capt})) ' m'])
    toto=[toto p{nn}];
    datetick('x')
end

linkaxes(toto, 'xy')

%% Celerity plot
for nn=1:2
    celerity{nn}=sndspd(sal{nn},temp{nn},mean(depth{nn}));
end

sal_moy=mean([sal{1} ; sal{2}]);
for nn=3:Nz+2
    celerity{nn}=sndspd(sal_moy,temp{nn},mean(depth{nn}));
end


toto=[];
figure
for nn=1:Nz+2
    ind_capt=b(nn);
    p{nn}=subplot(Nz+2,1,nn);
    plot(time_num{ind_capt},celerity{ind_capt})
    grid on
    title(['Mooring #' mooring ' - Depth = ' num2str(mean(depth{ind_capt})) ' m'])
    toto=[toto p{nn}];
    datetick('x')
end

linkaxes(toto, 'xy')

%% Temp profile

day=[2017,2,1];
day_num=datenum(day);

figure
for nn=1:Nz+2
   [~ , a]=min(abs(time_num{nn}-day_num));
   temp_z(nn)=temp{nn}(a);
   plot(temp_z(nn),mean(depth{nn}), '*k')
   hold on
end
xlabel('temperature (degree C)')
ylabel('depth (m)')
axis ij
title(['Mooring ' mooring ' ; ' datestr(day_num)])

grid on

