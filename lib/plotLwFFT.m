function [zMeterRel,zMeterUnrel] = plotLwFFT(ax, header, data, varargin)
% This assumes .lw6 files. No FFT or SNR subtraction is performed here, we
% assume data already containts spectra that are ready for plotting. 
hold(ax,'on'); 

% test if data is 1-dimensional 
if ~iscolumn((squeeze(data)))
    warning('data for FFT plotting is not 1-D...problems await!')
end

ci = []; 
if any(strcmpi(varargin,'ci')) 
    ci = varargin{find(strcmpi(varargin,'ci'))+1}; 
    ci = squeeze(ci); 
end

par = getParams(); 

% get frequencies of interest 
maxfreqidx = round(par.maxfreqlim/header.xstep)+1; 
frex_idx = round(par.frex/header.xstep)+1; 

freq = [0:header.xstep:par.maxfreqlim]; 

mX = squeeze(data); 
mX(1) = 0; 
mX = mX(1:maxfreqidx); 
if size(mX,1)>1 & size(data,2)==1
    mX = mX'; 
end

if ~isempty(ci)
    ci = ci(1:maxfreqidx); 
    if size(ci,1)>1 & size(data,2)==1
        ci = ci'; 
    end
   fill(ax, [freq,flip(freq)], [zeros(size(ci)),flip(mX+ci)], ...
        par.col_neutral, ...
        'facealpha',0.3, ...
        'LineStyle','none');
    
end

% calculate amplitudes and zscores 
amps = mX(frex_idx);
z = zscore(amps); 
zMeterRel = mean(z(par.idx_meterRel)); 
zMeterUnrel = mean(z(par.idx_meterUnrel)); 

stem(ax, freq,mX,'color',par.col_neutral,'marker','none','linew',par.linew_fft); 

stem(ax, freq(frex_idx(par.idx_meterRel)), mX(frex_idx(par.idx_meterRel)), ...
    'color',par.col_meterRel,...
    'marker','none',...
    'linew',par.linew_fft); 

stem(ax, freq(frex_idx(par.idx_meterUnrel)), mX(frex_idx(par.idx_meterUnrel)), ...
    'color',par.col_meterUnrel,...
    'marker','none',...
    'linew',par.linew_fft); 

%     yMax = round(1.1*max(mX)*par.prec)/par.prec; 
%     yMin = -0.05*yMax; 

ax = gca;
ax.Color = 'none'; 
ax.XTick = []; 
ax.XLim = [0,par.maxfreqlim]; 
%     ax.YLim = [yMin,yMax]; 
ax.YTick = [0,ax.YLim(2)]; 
ax.XAxis.Visible = 'off'; 
ax.FontSize = par.fontsize; 

