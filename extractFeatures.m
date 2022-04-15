function extractFeatures(roi)

par = getParams(); 

% we'll extract all harmonics up to 30Hz and process this then in R
frex_to_extract = par.frex_to_extract; 
idx_meterRel = [1:length(frex_to_extract)]; 
idx_meterUnrel = []; 

coch = load(fullfile(par.stim_path,'Slaney_128coch_meddis_fft_meanAcrossCF.mat')); 

fpath = par.features_path; 
fname = sprintf('features_roi-%s.csv',roi); 
fid = fopen(fullfile(fpath,fname),'w'); 
fprintf(fid,'subject,rhythm,tone,freq,isMeterRel,amp_eeg,amp_coch\n'); 

for iSub=1:length(par.subjects)

    for iRhythm=1:2

        for iTone=1:2

            subject = par.subjects(iSub); 
            tone = par.tones{iTone}; 
            rhythm = par.rhythms{iRhythm}; 

            %% EEG 
            
            fpath = fullfile(par.deriv_path, ...
                sprintf('roi-%s/sub-%03d',roi,subject)); 

            % load fft for the selected ROI 
            fname = sprintf('sub-%03d_rhythm-%s_tone-%s_snr-%d-%d_FFT',...
                subject,rhythm,tone,par.snr_bins_eeg(1), par.snr_bins_eeg(2)); 

            [header_fft,data_fft] = CLW_load(fullfile(fpath,fname)); 
            
            [~,~,z_eeg,~,~,amps_eeg] = getZ(squeeze(data_fft), frex_to_extract, ...
                idx_meterRel, idx_meterUnrel, 'xstep', header_fft.xstep); 

            %% coch

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

            [~,~,z_coch,~,~,amps_coch] = getZ(coch.res_all(idx,:), frex_to_extract, ...
                idx_meterRel, idx_meterUnrel, 'xstep', coch.freq(2)); 

            %% write to csv

            for fi=1:length(frex_to_extract)

                if ismember(fi,par.idx_meterRel)
                    isMeterRel = 1; 
                else
                    isMeterRel = 0; 
                end

                fprintf(fid,'%03d,%s,%s,%.3f,%d,%f,%f\n', ...
                    subject, rhythm, tone, frex_to_extract(fi), isMeterRel, ...
                    amps_eeg(fi), amps_coch(fi)); 

            end

        end
    end

end

fclose(fid); 






