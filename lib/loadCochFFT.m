function [header, data, freq] = loadCochFFT(rhythm, tone, varargin)

par = getParams(); 

fft_type = '';
if any(strcmpi(varargin, 'db'))
    fft_type = '_dB';
end

if any(strcmpi(varargin, 'db'))
    fft_type = '_dB'; 
    fprintf('loading cochlear model DECIBEL spectra...\n'); 
elseif any(strcmpi(varargin, 'pow'))
    fft_type = '_pow'; 
    fprintf('loading cochlear model POWER spectra...\n'); 
else
    fprintf('loading cochlear model AMPLITUDE spectra...\n'); 
    fft_type = ''; 
end
fname = sprintf('Slaney_128coch_meddis_fft%s_meanAcrossCF.mat', fft_type);
coch = load(fullfile(par.stim_path, fname));

maxfreqlim = 31;
if any(strcmpi(varargin, 'maxfreqlim'))
    maxfreqlim = varargin{find(strcmpi(varargin, 'maxfreqlim')) + 1};
end

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

row_idx = ~cellfun(@isempty, regexp(coch.cond_names, sprintf('^%s_',tone_code)), 'uni',1) & ...
          ~cellfun(@isempty, regexp(coch.cond_names, sprintf('_%s',rhythm_code)), 'uni',1);
row_idx = find(row_idx); 
assert(length(row_idx)==1)

maxfreqidx = dsearchn(coch.freq', maxfreqlim); 
mX = coch.mX(row_idx, 1:maxfreqidx); 
data = []; 
data(1,1,1,1,1,:) = mX; 

freq = coch.freq;
header = []; 
header.datasize = size(data); 
header.xstep = freq(2) - freq(1); 
