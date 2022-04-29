function [res, fs, t] = loadCochERP(rhythm,tone,varargin)

par = getParams(); 

if strcmpi(tone,'high')
    tone_code = 'H'; 
else
    tone_code = 'L'; 
end
if strcmpi(rhythm,'unsync')
    rhythm_code = 'unsyncopated'; 
else
    rhythm_code = 'syncopated'; 
end

coch = load(fullfile(par.stim_path,'Slaney_128coch_meddis_timeDomain_meanAcrossCF.mat')); 
fs = coch.fs; 

row_idx = ~cellfun(@isempty, regexp(coch.cond_names, sprintf('^%s_',tone_code)), 'uni',1) & ...
          ~cellfun(@isempty, regexp(coch.cond_names, sprintf('_%s',rhythm_code)), 'uni',1);

row_idx = find(row_idx); 
assert(length(row_idx)==1)


if any(strcmpi(varargin,'cycle'))
    [t, res] = getCycleErp(coch.hc(row_idx,:), coch.fs, 2.4); 
else
    res = coch.hc(row_idx,:);
    t = [0:length(res)-1]/coch.fs; 
end

res = res - min(res); 
res = res / max(res); 
% res = res - mean(res); 



