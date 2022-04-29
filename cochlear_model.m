
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
        t = [0:length(s)-1]/fs; 
        
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

        hN = length(hc(1,:))/2+1;
        freq = linspace(0, fs/2, hN);
        mX = abs(fft(hc,[],2)) / size(hc,2);
        mX = mX(:,1:hN);
        
        hc_all{end+1} = sum(hc,1);
        mX_all{end+1} = mean(mX,1);
        cond_names{end+1} = sprintf('%s_%s',tones{iTone}, rhythms{iRhythm});

    end
end

%% SAVE

hc = cat(1,hc_all{:});
mX = cat(1,mX_all{:});

fname = 'Slaney_128coch_meddis_fft_meanAcrossCF'; 
save(fullfile(par.stim_path,fname),'mX','cond_names','freq'); 

fname = 'Slaney_128coch_meddis_timeDomain_meanAcrossCF'; 
save(fullfile(par.stim_path,fname),'hc','cond_names','t','fs'); 










