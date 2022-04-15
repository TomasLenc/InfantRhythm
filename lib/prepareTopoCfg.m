function [cfg,d] = prepareTopoCfg(amps, chan_labels, lims, layout, varargin)

d = struct;
d.avg = amps; 
d.time = 1;
d.label = chan_labels;
d.dimord = 'chan_time';

cfg = [];                            
cfg.xlim = [1 1];   % time to plot             
cfg.zlim = lims;  % amplitude range              
cfg.layout = layout;     
cfg.marker = 'off'; 
cfg.comment = ' ';   
cfg.commentpos = 'title'; 
cfg.gridscale = 500;   
cfg.style = 'straight';

if any(strcmpi(varargin,'chan_sel')) 
    chan_sel = varargin{find(strcmpi(varargin,'chan_sel'))+1}; 
    cfg.highlight = 'on'; 
    cfg.highlightsymbol = '.'; 
    cfg.highlightsize = 5; 
    cfg.highlightcolor = 'k'; 
    cfg.highlightchannel = chan_sel; 
end

end