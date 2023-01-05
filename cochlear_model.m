
clear variables

par = getParams(); 

tones = {'L','H'}; 
rhythms = {'unsyncopated','syncopated'}; 

frex = 1/(12*0.2) * [1:12];
low_freq = 100;
nchan = 128;


%%
mX_all = {}; 
hc_all = {}; 
cond_names = {}; 

for iTone=1:length(tones)
    
    for iRhythm=1:length(rhythms)

        d = dir(fullfile(par.stim_path, ...
            sprintf('%s_%s_*.wav', tones{iTone}, rhythms{iRhythm})));
        assert(length(d)==1)
        disp(d.name); 
        
        [s,fs] = audioread(fullfile(d.folder, d.name));
        
        cfArray = ERBSpace(low_freq, fs/2, nchan);
        fcoefs = MakeERBFilters(fs, nchan, low_freq);
        coch = ERBFilterBank(s, fcoefs);

        % Meddis hair cell
        hc = MeddisHairCell(coch/max(max(coch))*80, fs, 1); 
        for j=1:size(coch,1)
            % low-pass filter
            c = hc(j,:);
            c=filter([1],[1 -.99],c);
            hc(j,:)=c;        
        end
        
        % butter filter to make it equivalent to the neural data!
        [b, a] = butter(...
            par.filter_order / 2,...
            [par.filter_low_cutoff, par.filter_high_cutoff] ./ (fs/2), ...
            'bandpass'...
            ); 
        
        for j=1:size(hc, 1)
            hc(j, :) = filtfilt(b, a, hc(j, :));
        end
        
        % cut off first and last few cycles to get rid of edge artifacts from
        % filter 
        assert(length(s)/fs == par.trial_duration)
        
        edge_rm_duration = par.cycle_duration * 2; 
        
        idx_start = round(edge_rm_duration * fs); 
        N = round((par.trial_duration - 2 * edge_rm_duration) * fs); 
        
        hc = hc(:, idx_start+1 : idx_start+N); 
        
        assert(size(hc, 2)/fs == par.trial_duration - 2 * edge_rm_duration)
        
        t = [0 : size(hc, 2)-1]/fs; 
        
%         % visual check for filtering artifacts
%         figure
%         plot(t, sum(hc, 1), 'linew', 1.5)
%         xlim([0, 20])
%         xlim([40, Inf])                
        
        % FFT
        hN = length(hc(1,:))/2+1;
        freq = linspace(0, fs/2, hN);
        mX = abs(fft(hc,[],2)) / size(hc,2);
        mX(1) = 0; 
        mX = mX(:,1:hN);
        
        hc_all{end+1} = sum(hc,1);
        mX_all{end+1} = mean(mX,1);
        cond_names{end+1} = sprintf('%s_%s',tones{iTone}, rhythms{iRhythm});

    end
end

%% SAVE

hc = cat(1, hc_all{:});
fname = 'Slaney_128coch_meddis_timeDomain_meanAcrossCF'; 
save(fullfile(par.stim_path,fname),'hc','cond_names','t','fs'); 

mX = cat(1, mX_all{:});
fname = 'Slaney_128coch_meddis_fft_meanAcrossCF'; 
save(fullfile(par.stim_path,fname),'mX','cond_names','freq'); 

mX = cat(1, mX_all{:}) .^ 2; 
fname = 'Slaney_128coch_meddis_fft_pow_meanAcrossCF'; 
save(fullfile(par.stim_path,fname),'mX','cond_names','freq'); 

mX = 20 * log10(cat(1, mX_all{:})); 
fname = 'Slaney_128coch_meddis_fft_dB_meanAcrossCF'; 
save(fullfile(par.stim_path,fname),'mX','cond_names','freq'); 










