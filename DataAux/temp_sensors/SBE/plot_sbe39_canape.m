% --------------------------------------------------
%   plot_sbe39_canape
%
%% --------------------------------------------------

save_plots = 1;

%% --------------------------------------------------
D = load('sbe3078.dat');

tstart = datenum(2011,4,8,4,0,0);
tend = datenum(2011,5,2,12,0,0);

T3078 = D(:,1);
P3078 = D(:,2);
day = D(:,3);
mon = D(:,4);
yr = D(:,5);
hr = D(:,6);
mm = D(:,7);
ss= D(:,8);

Tm3078 = datenum(yr,mon,day,hr,mm,ss);
Z3078 = sw_dpth(P3078,42);

figure(11);
plot(Tm3078,Z3078);

datetick; axis('ij');
grid on;
ylim([139 144]);
xlabel('Date 2016');
ylabel('Depth (m)');
title('Canape sbe39 Depth #3078');

if save_plots
    print -dpng sbe3078_z.png;
end

figure(12);
plot(Tm3078,T3078);
 datetick; 
grid on;
ylim([-2.5 0]);
xlabel('Date 2016');
ylabel('Temp (deg C)');
title('Canape sbe39 Temperature #3078');

if save_plots
    print -dpng sbe3078_t.png;
end
%% ----------------------------------------------------------
D = load('sbe3119.dat');

T3119 = D(:,1);
P3119 = D(:,2);
day = D(:,3);
mon = D(:,4);
yr = D(:,5);
hr = D(:,6);
mm = D(:,7);
ss= D(:,8);

Tm3119 = datenum(yr,mon,day,hr,mm,ss);
Z3119 = sw_dpth(P3119,42);

figure(13);
plot(Tm3119,Z3119);

datetick; axis('ij');
grid on;
ylim([138 145]);
xlabel('Date 2016');
ylabel('Depth (m)');
title('Canape sbe39 Depth #3119');

if save_plots
    print -dpng sbe3119_z.png;
end

figure(14);
plot(Tm3119,T3119);
 datetick; 
grid on;
ylim([-2.5 0]);
xlabel('Date 2016');
ylabel('Temp (deg C)');
title('Canape sbe39 Temperature #3119');

if save_plots
    print -dpng sbe3119_t.png;
end
%% ----------------------------------------------------------

D = load('sbe3123.dat');

T3123 = D(:,1);
P3123 = D(:,2);
day = D(:,3);
mon = D(:,4);
yr = D(:,5);
hr = D(:,6);
mm = D(:,7);
ss= D(:,8);

Tm3123 = datenum(yr,mon,day,hr,mm,ss);
Z3123 = sw_dpth(P3123,42);

figure(15);
plot(Tm3123,Z3123);

datetick; axis('ij');
grid on;
ylim([136 143]);
xlabel('Date 2016');
ylabel('Depth (m)');
title('Canape sbe39 Depth #3123');

if save_plots
    print -dpng sbe3123_z.png;
end

figure(16);
plot(Tm3123,T3123);
 datetick; 
grid on;
ylim([-2.5 0]);
xlabel('Date 2016');
ylabel('Temp (deg C)');
title('Canape sbe39 Temperature #3123');

if save_plots
    print -dpng sbe3123_t.png;
end
%% ----------------------------------------------------------

D = load('sbe3125.dat');

T3125 = D(:,1);
P3125 = D(:,2);
day = D(:,3);
mon = D(:,4);
yr = D(:,5);
hr = D(:,6);
mm = D(:,7);
ss= D(:,8);

Tm3125 = datenum(yr,mon,day,hr,mm,ss);
Z3125 = sw_dpth(P3125,42);

figure(15);
plot(Tm3125,Z3125);

datetick; axis('ij');
grid on;
ylim([137 144]);
xlabel('Date 2016');
ylabel('Depth (m)');
title('Canape sbe39 Depth #3125');

if save_plots
    print -dpng sbe3125_z.png;
end

figure(16);
plot(Tm3125,T3125);
 datetick; 
grid on;
ylim([-2.5 0]);
xlabel('Date 2016');
ylabel('Temp (deg C)');
title('Canape sbe39 Temperature #3125');

if save_plots
    print -dpng sbe3125_t.png;
end
%% ----------------------------------------------------------

D = load('sbe3130.dat');

T3130 = D(:,1);
P3130 = D(:,2);
day = D(:,3);
mon = D(:,4);
yr = D(:,5);
hr = D(:,6);
mm = D(:,7);
ss= D(:,8);

Tm3130 = datenum(yr,mon,day,hr,mm,ss);
Z3130 = sw_dpth(P3130,42);

figure(15);
plot(Tm3130,Z3130);

datetick; axis('ij');
grid on;
ylim([135 137]);
xlabel('Date 2016');
ylabel('Depth (m)');
title('Canape sbe39 Depth #3130');

if save_plots
    print -dpng sbe3130_z.png;
end

figure(16);
plot(Tm3130,T3130);
 datetick; 
grid on;
ylim([-2.5 0]);
xlabel('Date 2016');
ylabel('Temp (deg C)');
title('Canape sbe39 Temperature #3130');

if save_plots
    print -dpng sbe3130_t.png;
end
%% ----------------------------------------------------------



