close all
clear all
clc

close all
clear all
clc

%%%% wave
c=1440;
% f_=[350 600 900];
f_=350;
Nf=length(f_);

%%%% array
d=2.5;
n=4;

%%% look direction
t_s=90/180*pi;

t=linspace(-pi, pi, 250);

for ff=1:Nf
    f=f_(ff);    
    l=c/f;

    b=(sin(n*pi*d/l*sin(t)) ./ ( n*sin(pi*d/l*sin(t)) ));


    k=2*pi*f/c;
    b2=sin(n/2*d*k*( sin(t_s) - sin(t) ) ) ./ sin(1/2*d*k*( sin(t_s) - sin(t) ) );

%     figure(1)
%     plot(t/pi*180,10*log10(abs(b2).^2/max(abs(b2.^2))))
%     xlim([-90 90])
%     grid on
%     ylim([-45 5])
%     title('Broadside at 0 deg')
%     hold on

    % 
    figure
    polarplot(t,10*log10(abs(b2).^2/max(abs(b2.^2))), 'linewidth',2)
    rlim([-45 0])
    hold on
    title('Beampattern: up')
end


%%% look direction
t_s=0/180*pi;

t=linspace(-pi, pi, 250);

for ff=1:Nf
    f=f_(ff);    
    l=c/f;

    b=(sin(n*pi*d/l*sin(t)) ./ ( n*sin(pi*d/l*sin(t)) ));


    k=2*pi*f/c;
    b2=sin(n/2*d*k*( sin(t_s) - sin(t) ) ) ./ sin(1/2*d*k*( sin(t_s) - sin(t) ) );

%     figure(1)
%     plot(t/pi*180,10*log10(abs(b2).^2/max(abs(b2.^2))))
%     xlim([-90 90])
%     grid on
%     ylim([-45 5])
%     title('Broadside at 0 deg')
%     hold on

    % 
    figure
    polarplot(t,10*log10(abs(b2).^2/max(abs(b2.^2))), 'linewidth',2)
    rlim([-45 0])
    title('Beampattern: horizontal')
    hold on
end

% return
% 
% %%
% 
% sensorArrayAnalyzer;
% 
% 
% antenna = phased.ULA('NumElements',n,'ElementSpacing',d);
% 
% f_test=[f];
% azimuth=[-90:90];
% elevation=0;
% figure
% pattern(antenna,f_test,azimuth,elevation,'PropagationSpeed',c,...
%     'CoordinateSystem','rectangular', 'Type','powerdb');
% xlim([-90 90])
% ylim([-45 5])
% 
% figure
% pattern(antenna,f_test,azimuth,elevation,'PropagationSpeed',c,...
%     'CoordinateSystem','polar', 'Type','powerdb');
% 
