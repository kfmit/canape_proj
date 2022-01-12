addpath('/home/kfung/Downloads/CANAPE/mat_files/')
load('ANL_SHRU1.mat')
c = 1422;
f0 = 300;
fs = 3906.2;
p_ = [0 2.5 5 7.5];           % remember this is vertical
% checks
if ~iscolumn(p_)
    p_=transpose(p_);
end
n_chan = 4;

% Conventional Beamforming
% the general stuff: weighting, plotting, looking
lambda = c/f0;                              % wavelenght
theta=0:0.2:180;                            % angles to look on
Ntheta=length(theta);                       % Ntheta 
theta_t = 100;
v_t=exp(1i*cosd(theta_t)*p_*2*pi*f0/c);             % target at 95
    

for i_snap=1:5                       % is the nested forloop necessary
    s = SPL_raw(:,i_snap);                   % use one column of x_in as signal s that is 10x1
    
    for i=1:Ntheta
        k = -2*pi/lambda*cosd(theta(i));        % wavenumber
        v = exp(-1i*p_* k); % arraYy manifold vector according to k
        w_T = w_.*v;                         % weight times array manifold
        Y_CBF(i_snap,i) = w_T*s';                % element i is targ weight times vector
    end

end

% Plot conventional
figure(1)
plot(theta,abs(B_CBF))
% hold on
% plot(theta,abs(B_targ))