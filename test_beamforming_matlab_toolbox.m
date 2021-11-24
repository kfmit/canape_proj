clear all
close all
clc

fs=4e3;
c=1440;

%% Create array signal
t=0:1/fs:0.2;
N=length(t);

f0=300;
x=zeros(1,N);
h=hamming(floor(5/f0*fs)).';
deb=250;
x(deb:deb+length(h)-1)=h.*cos(2*pi*f0*t(1:length(h)));

noise=randn(1,N);

% xn=x+0.2*noise;
% plot(t,xn)


theta=45;
theta_look=45;

dz=2.5;
z=0:dz:7.5;
Nz=length(z);


xf=fft(x);
f=(0:N-1)*fs/N;
k=2*pi*f/c;

xf_vla=zeros(Nz,N);
for zz=1:Nz
    xf_vla(zz,:)=xf.*exp(-1j*k*z(zz)*cosd(theta));
end
x_vla=ifft(xf_vla,N,2, 'symmetric');   

% figure
% for zz=1:Nz
%     plot(t,x_vla(zz,:))
%     hold on
%     grid on
% end 
% legend
% return

noise=0.5*randn(size(x_vla));
x_vla_n=x_vla+noise;

figure
for zz=1:Nz
    subplot(121)
    plot(t,x_vla(zz,:))
    title('Signal on array, no noise')
    hold on
    subplot(122)
    plot(t,x_vla_n(zz,:))
    title('Signal on array, with noise')
    hold on
end


%% Beamform



sig=x_vla_n;

tic
sig_f=fft(sig,N,2);
d=zeros(size(sig_f));
for zz=1:Nz
    d(zz,:)=exp(1j*k*z(zz)*cosd(theta_look));
end
tic
sig_f_beam=sig_f.*d;
sig_t_beam=ifft(sig_f_beam,N,2, 'symmetric');
sig_ok=sum(sig_t_beam,1)/Nz;
t_ju=toc




%% Use Matlab toolbox

antenna = phased.ULA('NumElements',Nz,'ElementSpacing',dz);

f_test=[250:50:350];
azimuth=[-180:180];
elevation=0;
figure
pattern(antenna,f_test,azimuth,elevation,'PropagationSpeed',c,...
    'CoordinateSystem','rectangular', 'Type','directivity');


% [PAT,AZ_ANG,EL_ANG] = pattern(antenna,f_test,azimuth,elevation,'PropagationSpeed',c,...
%     'CoordinateSystem','rectangular', 'PlotStyle', 'waterfall', 'Type','directivity');

figure,
imagesc(f_test,azimuth,PAT)
xlabel('Frequency (Hz)')
ylabel('Angle (degree)')


collector = phased.WidebandCollector('Sensor',antenna,'PropagationSpeed',c,'SampleRate',fs,...
    'NumSubbands',N,'ModulatedInput', false);
%%%%% matlab has broadside at 0 degree
incidentAngle=[theta-90 ; 0];

x_col=x.';
x_array = collector(x_col,incidentAngle);

figure
p1=subplot(211);
plot(t,x_array)
title('Signal on array, manual')
p2=subplot(212);
plot(t,x_vla.')
title('Signal on array, matlab')
linkaxes([p1  p2], 'xy')


x_array_n=x_array+noise.';

%%% Matlab has broadside at 0 degree
beamform_angle=[theta_look-90 ; 0];

beamformer = phased.TimeDelayBeamformer('SensorArray',antenna, ...
    'Direction',beamform_angle, 'PropagationSpeed',c,'SampleRate',fs, 'WeightsOutputPort',true);


tic
[y,w] = beamformer(x_array_n);
t_mat=toc

figure
subplot(211)
plot(t,sig_ok,t,real(y))
grid on
title('recovered signal')
legend('manual', 'matlab')
subplot(212)
plot(t,x_vla_n(1,:))
grid on
title('original signal')

