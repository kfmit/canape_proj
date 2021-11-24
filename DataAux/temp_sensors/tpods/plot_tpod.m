function plot_tpod(fname)
% -----------------------------------------------
%  plot_tpod
%
%  Must be prepared for matlab input
%    meta data commented with %
%    no /'s in date, should be space
%    no :'s in time, should be space
%   real data must have .   not  ,
% -------------------------------------------------

if nargin < 1
    error('   USAGE: plot_tpod(FILENAME)');
end

% load tpod file
tdata = load(fname);



day = tdata(:,2);
mon =  tdata(:,3);
yr =  2000+tdata(:,4);
hr = tdata(:,5);
mins = tdata(:,6);
sec = tdata(:,7);
temp = tdata(:,8);

tm = datenum(yr,mon,day,hr,mins,sec);

figure(1);
plot(tm,temp); grid on
% set(gca,'xlim',[tmin tmax]);
% set(gca,'ylim',[16 17]);
datetick('x','keeplimits');

title([fname]);
xlabel('Date');
ylabel('Temperature (deg C)')









