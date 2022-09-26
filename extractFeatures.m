function extractFeatures(roi)

par = getParams(); 

% we'll extract all harmonics up to 30Hz and process this then in R
frex_to_extract = par.frex_to_extract; 
idx_meterRel = [1:length(frex_to_extract)]; 
idx_meterUnrel = []; 

fpath = par.features_path; 
fname = sprintf('features_roi-%s.csv', roi); 
fid = fopen(fullfile(fpath,fname),'w'); 
fprintf(fid,'subject,rhythm,tone,freq,isMeterRel,amp_eeg,db_eeg,amp_coch,pow_coch,db_coch\n'); 

for iSub=1:length(par.subjects)
    
    subject = par.subjects(iSub); 
    fprintf('\n\n----------------\nsub-%d\n', subject); 

    for iRhythm=1:2

        rhythm = par.rhythms{iRhythm}; 
        
        for iTone=1:2

            tone = par.tones{iTone}; 

            %% EEG 
            
            fpath = fullfile(par.deriv_path, ...
                sprintf('roi-%s/sub-%03d',roi,subject)); 

            % load fft magnitudes for the selected ROI 
            fname = sprintf('sub-%03d_rhythm-%s_tone-%s_snr-%d-%d_FFT',...
                subject,rhythm,tone,par.snr_bins_eeg(1), par.snr_bins_eeg(2)); 

            [header_fft, data_fft] = CLW_load(fullfile(fpath,fname)); 
            
            [~,~,~,~,~, amps_eeg] = getZ(...
                squeeze(data_fft), frex_to_extract, ...
                idx_meterRel, idx_meterUnrel,...
                'xstep', header_fft.xstep...
                ); 
            
            % load fft decibels for the selected ROI 
            fname = sprintf('sub-%03d_rhythm-%s_tone-%s_snr-%d-%d_FFT_dB',...
                subject,rhythm,tone,par.snr_bins_eeg(1), par.snr_bins_eeg(2)); 

            [header_fft, data_fft] = CLW_load(fullfile(fpath,fname)); 
            
            [~,~,~,~,~, db_eeg] = getZ(...
                squeeze(data_fft), frex_to_extract, ...
                idx_meterRel, idx_meterUnrel,...
                'xstep', header_fft.xstep...
                ); 

            %% coch

            [header_coch, data_coch] = loadCochFFT(rhythm, tone);

            [~,~,~,~,~, amps_coch] = getZ(...
                squeeze(data_coch), frex_to_extract, ...
                idx_meterRel, idx_meterUnrel, ...
                'xstep', header_coch.xstep); 

            [header_coch, data_coch] = loadCochFFT(rhythm, tone, 'pow');

            [~,~,~,~,~, pow_coch] = getZ(...
                squeeze(data_coch), frex_to_extract, ...
                idx_meterRel, idx_meterUnrel, ...
                'xstep', header_coch.xstep); 

            [header_coch, data_coch] = loadCochFFT(rhythm, tone, 'db');

            [~,~,~,~,~, db_coch] = getZ(...
                squeeze(data_coch), frex_to_extract, ...
                idx_meterRel, idx_meterUnrel, ...
                'xstep', header_coch.xstep); 
            
            
            %% write to csv

            for fi=1:length(frex_to_extract)

                if ismember(fi,par.idx_meterRel)
                    isMeterRel = 1; 
                else
                    isMeterRel = 0; 
                end

                fprintf(...
                    fid,'%03d,%s,%s,%.3f,%d,%f,%f,%f,%f,%f\n', ...
                    subject, rhythm, tone, frex_to_extract(fi), isMeterRel, ...
                    amps_eeg(fi), ...
                    db_eeg(fi), ...
                    amps_coch(fi), ...
                    pow_coch(fi), ...
                    db_coch(fi) ...
                    ); 

            end

        end
    end

end

fclose(fid); 






