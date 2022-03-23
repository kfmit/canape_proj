function [alpha1,alpha2,alpha3] = attenuation(f,c,type)
%UNTITLED6 Summary of this function goes here
%   specify unit as 1 = dB/m, 2 = dB/lambda
% c is the sound speed vector
lambda = c./f;

% conditions where certain equations may be used
alpha_f0 = 0.26; % if f0 is 1000
f0=1000;
u = 1.8; % if f<1kHz

if f<=1000
    disp('f<1kHz, all 3 alphas')
    alpha1 = 0.37*(f/1000)^1.80;
    alpha2 = 0.33*(f/1000)^1.86;
    alpha3 = alpha_f0*(f/f0)^u;
elseif f>1000
    f_kHz=f/1000;
    alpha2 = 0.33*f_kHz^1.86;
    alpha1 = nan;
    alpha3 = nan;
    disp('No alpha 1 or alpha3')
end

if type==1
    disp('Units in dB/m')
elseif type==2
    disp('Units in dB/lambda')
    alpha1 = lambda.*alpha1;
    alpha2 = lambda.*alpha2;
    alpha3 = lambda.*alpha3;
end
% end of function 
end