%%%% Percentile
perc_plot=[1 5 50 95 99];

%%%% dB axis
mindB=min(min(A));
maxdB=max(max(A));
% mindB=db_lim_plot(1);
% maxdB=db_lim_plot(2);
hind = 0.1;               % histogram bin width for probability densities (PD)
dbvec = mindB:hind:maxdB; % dB values at which to calculate empirical PD

%%% Compute probability densities
d = hist(A,dbvec)/(hind*size(A,1));  %SPD array
d(d == 0) = NaN;   %suppress plotting of empty hist bins

%%% Compute percentiles
p=prctile(A,perc_plot,1);

%%% Compute RMS level
RMSlevel = 10*log10(mean(10.^(A/10)));

%%% Plot
g = pcolor(fPSD,dbvec,d); 
set(g,'LineStyle','none')
% set(gca, 'XScale', 'log')
colorbar
grid on
hold on
plot(fPSD, p, 'k', 'linewidth', 2)
hold on
plot(fPSD, RMSlevel, 'r', 'linewidth', 2)
xlabel('Frequency (Hz)')
ylabel('PSD (dB re 1 \muPa^2/Hz)')
ylabel(colorbar,'Empirical Probability Density','fontsize',14,'fontname','Arial')
% caxis([0 0.05])

if ~isempty(db_lim_plot)
    ylim(db_lim_plot)
end