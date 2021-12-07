%% Geographic Location of SHRUs

clear
clc
close all

% LATITUDES
shrulat = [54.4123 45.2347 40.6924 36.6582 45.4580]/60;         % N
shrulat = shrulat+72;                                           % latitude
shrulon =[159.018066666666667 158.2721 157.9108 157.5375 157.4874 ];    % W
shrudepth = [163 163 163 163 163];  % start depth of each thing
shrud_d = [0 2.5 5 7.5];
% water depth is NOT shruepth
w_depth = [302 308 303 301 445]*-1;                           % in m, doesn't correlate with paper??
shrutime = [21:04:00 23:37:00 02:04:00 04:21:00 23:28:00];
shrudate = ['24-Oct-2016' '24-Oct-2016' '25-Oct-2016' '25-Oct-2016' '21-Oct-2016'];
% Secondary structure
location = zeros(length(shrulat),3); % matrix of latitude, longitude
location(:,1) = shrulat;
location(:,2) = shrulon;
location(:,3) = w_depth;

% Tertiary structure
SHRUinfo = struct('lat',shrulat,'lon',shrulon,'depth', shrudepth,'w_depth', ...
                    w_depth,'time',shrutime,'date',shrudate,'space',shrud_d);

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

%% visualization of one nugget

figure
plot(SHRUinfo.lat(1),-1*SHRUinfo.space)

%% Resolution of Whole Array

figure(2)
hold on
scatter(SHRUinfo.lat(1),SHRUinfo.lon(1),'^','filled')
scatter(SHRUinfo.lat(2),SHRUinfo.lon(2),'^','filled')
scatter(SHRUinfo.lat(3),SHRUinfo.lon(3),'^','filled')
scatter(SHRUinfo.lat(4),SHRUinfo.lon(4),'^','filled')
scatter(SHRUinfo.lat(5),SHRUinfo.lon(5),'^','filled')
hold off

legend('SHRU1','SHRU2', 'SHRU3', 'SHRU4', 'SHRU5')
zlim([-500 0])
xlabel('Latitude ^\circ N')
ylabel('Longitude ^\circ W')
title('Hydrophone Location and Depth')

%% Beampattern of the WHOLE array


% in the horizontal direction: uses latitude


% in the veritcal direction: uses longitude (makes most sense)


% the whole 3d array: use the distance formula

%% Beamforming part of the array
% aka just one shru: sensors ar 2.5 m apaprt

n_chan = 4;     % in 
c = 1442;       % in channel speed
d = 2.5;        % in m
fs = 3906.2;    % sample freq
f0 = 300;       % Hz, what we're looking for 
p_ = [0 2.5 5 7.5];           % remember this is vertical
% checks
if ~iscolumn(p_)
    p_=transpose(p_);
end


% Conventional Beamforming
% the general stuff: weighting, plotting, looking
lambda = c/f0;                              % wavelenght
theta=0:0.2:180;                            % angles to look on
Ntheta=length(theta);                       % Ntheta 
theta_t = 100;
v_t=exp(1i*cosd(theta_t)*p_*2*pi*f0/c);             % target at 95
%p_ = % element locations established above

% Make a beampattern
B_CBF=zeros(1, Ntheta);                     % allocation
B_targ=zeros(1, Ntheta);                    % allocations
w_ = ones(n_chan,1);                      % weights for conventional
w_targ = w_.*v_t;                           % target weighting
    

for i=1:Ntheta
    k = -2*pi/lambda*cosd(theta(i));        % wavenumber
    v = exp(-1i*p_* k);                     % array manifold vector according to k
    B_CBF(i) = (w_'*v);                     % element i is weight times vector
    B_targ(i) = (w_targ'*v);                % element i is targ weight times vector
    %     Y_CBF(i_snap,i) = sum(w_'*s);
end

clear v k

% Plot conventional
figure(1)
plot(theta,abs(B_CBF))
% hold on
% plot(theta,abs(B_targ))
grid on
xlim([theta(1) theta(end)]);
ylabel('|\Upsilon| Output')
xlabel('\theta')
% legend('Broadside','Target')
title('Vertical: Conventional Beampattern |\Upsilon|')

figure
plot(theta,20*log(abs(B_CBF)))
ylabel('|\Upsilon| Output')
xlabel('\theta')
title('Vertical: Conventional Beampattern Power of |\Upsilon|')

%% Verical Array: Using the matlab toolbox
H = phased.ULA('NumElements',4,'ElementSpacing',2.5);%'ArrayAxis','z');
p_matlab =getElementPosition(H);
% it doesn't like when I set it to the z axis
viewArray(H)

f_test = [250 275 300 325 350];
figure
pattern(H,f_test,-180:180,0,'PropagationSpeed',c,'CoordinateSystem','rectangular','Type','powerdb','Normalize',true)

figure
plotResponse(H,f_test,c)


%% test beamforming on Horizontal

% find the distance between each hydrophone, make that the new p and new
% weight
% Handle as line array wher 5 is excluded?
dist_shrus = zeros(5,1);
for i=2:5
    degdist = distance('gc',SHRUinfo.lat(i),SHRUinfo.lon(i),SHRUinfo.lat(i-1),SHRUinfo.lon(i-1));
    dist_shrus(i) = deg2km(degdist);
end

%%% NEW p_

for ii=1:5
    p_(ii)=sum(dist_shrus(1:ii));
end

figure
scatter(p_,zeros(1,length(p_)),70, 'filled')
title('Geometry of Horizontal Array')
xlabel('meters')
%%
n_chan = 5;     % in 
c = 1442;       % in channel speed
d = 2.5;        % in m
fs = 3906.2;    % sample freq
%f0 = 300;       % Hz, what we're looking for 

% checks
if ~iscolumn(p_)
    p_=transpose(p_);
end

for n=1:length(f_test)
    f0=f_test(n);
% Conventional Beamforming
% the general stuff: weighting, plotting, looking
lambda = c/f0;                              % wavelenght
theta=0:0.2:180;                            % angles to look on
Ntheta=length(theta);                       % Ntheta 
theta_t = 100;
v_t=exp(1i*cosd(theta_t)*p_*2*pi*f0/c);             % target at 95
%p_ = % element locations established above

% Make a beampattern
B_CBF=zeros(1, Ntheta);                     % allocation
B_targ=zeros(1, Ntheta);                    % allocations
w_ = ones(n_chan,1);                      % weights for conventional
w_targ = w_.*v_t;                           % target weighting
    

for i=1:Ntheta
    k = -2*pi/lambda*cosd(theta(i));        % wavenumber
    v = exp(-1i*p_* k);                     % array manifold vector according to k
    B_CBF(i) = (w_'*v);                     % element i is weight times vector
    B_targ(i) = (w_targ'*v);                % element i is targ weight times vector
    %     Y_CBF(i_snap,i) = sum(w_'*s);
end

clear v k

% Plot conventional
figure(1)
plot(theta,abs(B_CBF))
hold on
% plot(theta,abs(B_targ))
grid on
end 
xlim([theta(1) theta(end)]);
ylabel('|\Upsilon| Output')
xlabel('\theta')
% legend('Broadside','Target')
title('Horizontal: Conventional Beampattern |\Upsilon|')
legend('250',' 275','300','325','350')

figure
plot(theta,20*log(abs(B_CBF)))
ylabel('|\Upsilon| Output')
xlabel('\theta')
title('Horizontal: Conventional Beampattern Power of |\Upsilon|')

%% Horizontal Array: Using the matlab toolbox

H = phased.ULA('NumElements',4,'ElementSpacing',p_);%'ArrayAxis','z');

% Handle as irregular array?

% handle as triangular array
