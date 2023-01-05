function [t, erp_avg, erp_sem, erp_all_foi] = getCycleErp(data, fs, cycle_duration, varargin)
% this function cuts data over the LAST dimension and averages the chunks 
% the input can be either 1D array, or letswave 6D data array

parser = inputParser; 

addParameter(parser, 'start_time', 0); 

parse(parser, varargin{:}); 

start_time = parser.Results.start_time; 


% fix 1D data to be a row vector (remember: time must be the last dimension!)
if size(data,1)>1 & size(data,2)==1
    data = data'; 
end

trial_duration = length(data)/fs; 

chunk_dur    = cycle_duration; 
chunk_n      = round(chunk_dur*fs); 
n_chunks     = floor((trial_duration - start_time) / chunk_dur - 1); 
n_foi        = size(data,5); 

erp = nan(n_chunks,n_foi,chunk_n); 

for i=1:n_chunks
    start_idx = round((start_time+(i-1)*chunk_dur)*fs);
    % go over last dimension (this makes it generalize across inputs of any
    % dimension)
    v = repmat({':'},ndims(data),1); 
    v{end} = start_idx+1:start_idx+chunk_n; 
    erp(i,:,:) = data(v{:}); 
end

% average across chunks 
erp_avg = squeeze(mean(erp,1))'; 
erp_sem = squeeze(std(erp,[],1)/sqrt(size(erp,1)))'; 

% if we have multiple frequency lines, save them in a separate field 
if all(size(erp_avg)>1)
    erp_all_foi = erp_avg; 
    erp_avg = []; 
else 
    erp_all_foi = []; 
end

t = [0:chunk_n-1] / fs; 
