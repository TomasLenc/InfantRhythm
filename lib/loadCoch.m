function [mX,xstep] = loadCoch(rhythm,tone)

par = getParams(); 

coch = load(fullfile(par.stim_path,'Slaney_128coch_meddis_fft_meanAcrossCF.mat')); 

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
fname = sprintf('%s_%s',tone_code,rhythm_code);
idx = find(strcmpi(fname,coch.rowNames)); 

mX = coch.res_all(idx,:); 
xstep = coch.freq(2)-coch.freq(1); 