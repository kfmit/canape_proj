% -----------------------------------------
%   plot_canape_tpods
%
%
% ----------------------------------------

save_plots = 0;

tmin = datenum(2016,10,19,0,0,0);
tmax = datenum(2017,10,19,0,0,0);

all_files = dir('*.out');

for ifile = 1:length(all_files)
    
    fname = all_files(ifile).name;
    plot_tpod(fname);
    
    S = sprintf('Canape2016 DRDC tpods %s',fname(1:end-4) );
    title(S)
    xlabel('Date (2016-2017)');
    ylabel('Temperature (deg C)');
    
    xlim([tmin tmax]);
    ylim([-2 2]);
    
       drawnow
       
    if save_plots
        S = sprintf('%s_temp.png',fname(1:end-4));
        print('-dpng', S);
    end
    
  
end
