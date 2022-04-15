function plotAmpAround(ax, amps, z, snr_bins, varargin)

prec = 2; 
if any(strcmp(varargin,'prec'))
    prec = varargin{find(strcmp(varargin,'prec'))+1}; 
end

yLimMin = floor(min(amps)*prec)/prec; 
yLimMax = ceil(max(amps)*prec)/prec; 
yLims = [yLimMin*0.8, yLimMax*1.2]; 

plot(ax,amps,...
    'color',[0.5,0.5,0.5], ...
    'linew',1.5); 

hold(ax,'on')

tmp = plot(ax,[snr_bins(2)+1,snr_bins(2)+1], yLims, ':r', 'linew',1.5); 
tmp.Color(4) = 0.5; 
txt = text(ax, snr_bins(2)*2, yLims(2), sprintf('z=%.1f',z)); 
txt.HorizontalAlignment = 'right'; 

ax.YLim = yLims; 
ax.Color = 'none'; 
ax.Visible = 'off';
