function [s,fs,t] = loadStim(rhythm,tone,varargin)

par = getParams(); 


if strcmpi(rhythm,'unsync')
    rhythm_code = 'unsyncopated'; 
else
    rhythm_code = 'syncopated'; 
end
if strcmpi(tone,'high')
    tone = 'H'; 
else
    tone = 'L'; 
end


fname = sprintf('%s_%s_*.wav',tone,rhythm_code); 
d = dir(fullfile(par.stim_path, fname)); 

[s,fs] = audioread(fullfile(d.folder,d.name)); 

if any(strcmpi(varargin,'cycle'))
   idx_start = round(par.erp_t_start_plot(rhythm) * fs) + 1; 
   idx_end = round(par.cycle_duration * fs); 
   s = s(idx_start : idx_end); 
end

t = [0:length(s)-1]/fs; 