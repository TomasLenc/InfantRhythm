function [zMeterRel,zMeterUnrel] = plotCochFFT(ax, x, fs, snr_bins)
% Plot FFT magnitude spectra, or autocorelation. 

par = getParams(); 

[mX,freq,amps,z,zMeterRel,zMeterUnrel] = getFFT( ...
    x, fs, par.frex, par.idx_meterRel, ...
    'maxfreqlim',par.maxfreqlim, 'snr_bins',snr_bins); 

frex_idx = round(par.frex/fs*length(x))+1; 

stem(freq,mX,'color',par.col_neutral,'marker','none','linew',par.linew_fft); 

hold on

stem(freq(frex_idx(par.idx_meterRel)), mX(frex_idx(par.idx_meterRel)), ...
    'color',par.col_meterRel,...
    'marker','none',...
    'linew',par.linew_fft); 

stem(freq(frex_idx(par.idx_meterUnrel)), mX(frex_idx(par.idx_meterUnrel)), ...
    'color',par.col_meterUnrel,...
    'marker','none',...
    'linew',par.linew_fft); 

ax = gca;
ax.Color = 'none'; 
ax.XTick = []; 
ax.XLim = [0,par.maxfreqlim]; 
ax.YTick = []; 
ax.XAxis.Visible = 'off'; 
ax.FontSize = par.fontsize; 

    