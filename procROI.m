function procROI(subject, roi)

par = getParams(); 

fprintf('sub-%03d  roi-%s\n',subject,roi); 

for iRhythm=1:2

    for iTone=1:2

        tone = par.tones{iTone}; 
        rhythm = par.rhythms{iRhythm}; 

        % load data
        fname = sprintf('AB %s %s P%03d*.lw6',tone,rhythm,subject); 
        d = dir(fullfile(par.preproc_path, fname)); 
        [header,data] = CLW_load(fullfile(d.folder,d.name)); 

        % rereference to mastoid (E57, E100)
        [header,data] = RLW_rereference(header,data,...
            'apply_list',{header.chanlocs.labels},'reference_list',par.ref_chans); 

        % filter for ERPs 
        [header_lp,data_lp] = RLW_butterworth_filter(header, data, ...
            'filter_type', 'lowpass', 'high_cutoff', 20, 'filter_order', 4); 
        
        % cut trials into chunks before averaging ? 
        if par.preproc_chunk_do
            [header_ep,data_ep,msg] = RLW_segmentation_chunk(header, data,...
                'chunk_onset',par.preproc_chunk_onset, ...
                'chunk_duration',par.preproc_chunk_dur, ...
                'chunk_interval',par.preproc_chunk_dur); 
            disp(msg); 
        else
            header_ep = header; 
            data_ep = data; 
        end

        % average epochs 
        [header_ep,data_ep] = RLW_average_epochs(header_ep, data_ep); 
        [header_lp,data_lp] = RLW_average_epochs(header_lp, data_lp); 

        % FFT
        [header_fft_noSNR, data_fft_noSNR] = RLW_FFT(header_ep, data_ep);     

        % SNR (2-5)
        [header_fft, data_fft] = RLW_SNR(header_fft_noSNR,data_fft_noSNR, ...
            'xstart', par.snr_bins_eeg(1), 'xend',par.snr_bins_eeg(2)); 
        
        % convert spectrum to dB
        header_fft_db = header_fft_noSNR;
        data_fft_db = 20 * log10(data_fft_noSNR); 
        
        % perform SNR noise subtraction 
        [header_fft_db, data_fft_db] = RLW_SNR(header_fft_db, data_fft_db, ...
            'xstart', par.snr_bins_eeg(1), 'xend',par.snr_bins_eeg(2)); 
        
        
        %% select channels 

        if strcmp(roi,'front')

            % pool channels (spectra averaged across 28 frontocentral channels (14 in each
            % hemisphere))
            chan_sel = par.front_chan'; 

        elseif strcmp(roi,'best')

            % laod data for selected channel 
            fpath = fullfile(par.deriv_path, ...
                sprintf('roi-best/sub-%03d',subject)); 
            fname = sprintf('chan_names.csv'); 
            chan_sel = readtable(fullfile(fpath,fname)); 
            chan_sel = chan_sel.label; 

        end
        
        % average roi channels in the time-domain (for ERPs)
        [header_time_roi, data_time_roi] = RLW_pool_channels(...
            header_lp, data_lp, chan_sel,...
            'mixed_channel_label', 'avg', 'keep_original_channels', false...
            ); 

        % average spectra (magnitudes) across roi channels 
        [header_fft_roi, data_fft_roi] = RLW_pool_channels(...
            header_fft, data_fft, chan_sel,...
            'mixed_channel_label', 'avg', 'keep_original_channels', false...
            ); 

        % average spectra (magnitudes) across roi channels (FFT without noise
        % subtraction)
        [header_fft_noSNR_roi, data_fft_noSNR_roi] = RLW_pool_channels(...
            header_fft_noSNR, data_fft_noSNR, chan_sel,...
            'mixed_channel_label','avg','keep_original_channels',false...
            ); 

        % average spectra (magnitudes) across roi channels (dB-transformed)
        [header_fft_db_roi, data_fft_db_roi] = RLW_pool_channels(...
            header_fft_db, data_fft_db, chan_sel,...
            'mixed_channel_label', 'avg', 'keep_original_channels', false...
            ); 
        
        
        %% save 
        
        fpath = fullfile(par.deriv_path, ...
            sprintf('roi-%s/sub-%03d',roi,subject)); 

        if ~isdir(fpath); mkdir(fpath); end 

        % if applicable, save (or re-save) selected channel names 
        if ~isempty(chan_sel)
            fname = sprintf('chan_names.csv'); 
            chan_table = cell2table(chan_sel, 'VariableNames',{'label'}); 
            writetable(chan_table,fullfile(fpath,fname)); 
        end    

        % time-domain 
        fname = sprintf('sub-%03d_rhythm-%s_tone-%s_time',subject,rhythm,tone); 
        header_time_roi.name = fname; 
        CLW_save(fpath, header_time_roi, data_time_roi); 

        % FFT no-SNR
        fname = sprintf('sub-%03d_rhythm-%s_tone-%s_snr-0-0_FFT',...
                            subject,rhythm,tone); 
        header_fft_noSNR_roi.name = fname; 
        CLW_save(fpath, header_fft_noSNR_roi, data_fft_noSNR_roi); 

        % FFT 
        fname = sprintf('sub-%03d_rhythm-%s_tone-%s_snr-%d-%d_FFT',...
            subject,rhythm,tone,par.snr_bins_eeg(1), par.snr_bins_eeg(2)); 
        header_fft_roi.name = fname; 
        CLW_save(fpath, header_fft_roi, data_fft_roi); 

        % FFT log-transformed
        fname = sprintf('sub-%03d_rhythm-%s_tone-%s_snr-%d-%d_FFT_dB',...
            subject,rhythm,tone,par.snr_bins_eeg(1), par.snr_bins_eeg(2)); 
        header_fft_db_roi.name = fname; 
        CLW_save(fpath, header_fft_db_roi, data_fft_db_roi); 
        
    end
end




