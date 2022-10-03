function procAllChan(subject)
% This function calculates FFT on preprocessed data across all channels. 
% Also, it calculates metrics across channels that will be needed to extract
% ROIs. 


par = getParams(); 

% load table with bad subjects 
fname = fullfile(par.deriv_path,'bad_subjects.csv'); 
if ~exist(fname)
    bad_subjects_table = cell2table(cell(0,1), 'VariableNames',{'subject'}); 
else
    bad_subjects_table = readtable(fname); 
end
    
fprintf('sub-%03d\n',subject); 

% allocate overall response magnitude - this will be used to find best
% channels across all conditions 
amps_sum = nan(par.n_chan,2,2); 

%% 

for iRhythm=1:2

    for iTone=1:2

        tone = par.tones{iTone}; 
        rhythm = par.rhythms{iRhythm}; 

        % load data
        fname = sprintf('AB %s %s P%03d*.lw6',tone,rhythm,subject); 
        d = dir(fullfile(par.preproc_path, fname)); 
        [header,data] = CLW_load(fullfile(d.folder,d.name)); 

        % rereference to mastoid (E57, E100)
        [header,data] = RLW_rereference(...
            header,data,...
            'apply_list',{header.chanlocs.labels}, ...
            'reference_list',par.ref_chans...
            ); 

        % cut trials into chunks before averaging ? 
        if par.preproc_chunk_do
            [header,data,msg] = RLW_segmentation_chunk(...
                header,data,...
                'chunk_onset',par.preproc_chunk_onset, ...
                'chunk_duration',par.preproc_chunk_dur, ...
                'chunk_interval',par.preproc_chunk_dur...
                ); 
            disp(msg); 
        end
        
        % average epochs 
        [header,data] = RLW_average_epochs(header,data); 

        fs = 1/header.xstep; 
        N = header.datasize(end); 

        % FFT
        [header_fft_noSNR, data_fft_noSNR] = RLW_FFT(header,data);             

        % SNR (2-5)
        [header_fft, data_fft] = RLW_SNR(header_fft_noSNR,data_fft_noSNR, ...
            'xstart', par.snr_bins_eeg(1), 'xend',par.snr_bins_eeg(2)); 

        % convert spectrum to dB
        header_fft_db = header_fft_noSNR;
        data_fft_db = 20 * log10(data_fft_noSNR); 
        
        % perform SNR noise subtraction 
        [header_fft_db, data_fft_db] = RLW_SNR(header_fft_db, data_fft_db, ...
            'xstart', par.snr_bins_eeg(1), 'xend',par.snr_bins_eeg(2)); 
        
        % infants showing large residual artifacts (amplitudes larger than 10 microV in 
        % more than 10 channels within one of the conditions) within the frequency range 
        % of interest (0.5 to 5 Hz) were then excluded from further analysis (5 infants over 20)
        idx_min = round(0.5/header_fft.xstep)+1; 
        idx_max = round(5/header_fft.xstep)+1; 
        d = any(squeeze(data_fft(:,:,1,1,1,idx_min:idx_max))>10, 2); 
        if sum(d)>10
            warning('too much noise in: sub-%03d %s %s',subject,rhythm,tone);
            if ~ismember(subject, bad_subjects_table.subject)
                bad_subjects_table.subject(end+1) = subject; 
            end
        else
            if ismember(subject, bad_subjects_table.subject)
                idx = find(bad_subjects_table.subject == subject); 
                bad_subjects_table = ...
                    bad_subjects_table(bad_subjects_table.subject ~= subject, :); 
            end
        end

        % get amplitudes for best-channel ROI selection 
        [~,~,~,amps,~] = getZ(...
            squeeze(data_fft), par.frex, ...
            par.idx_meterRel, par.idx_meterUnrel, ...
            'xstep', header_fft.xstep...
            ); 

        amps_sum(:,iRhythm,iTone) = sum(amps,2); 
        
        
        % save data for all channels
        % --------------------------

        fpath = fullfile(par.deriv_path, ...
            sprintf('roi-all/sub-%03d',subject)); 

        if ~isdir(fpath); mkdir(fpath); end 

        % without SNR 
        fname = sprintf('sub-%03d_rhythm-%s_tone-%s_snr-0-0_FFT',...
            subject, rhythm, tone); 
        header_fft_noSNR.name = fname; 
        CLW_save(fpath, header_fft_noSNR, data_fft); 

        % with SNR 
        fname = sprintf('sub-%03d_rhythm-%s_tone-%s_snr-%d-%d_FFT',...
            subject, rhythm, tone, par.snr_bins_eeg(1), par.snr_bins_eeg(2)); 
        header_fft.name = fname; 
        CLW_save(fpath, header_fft, data_fft); 

        % log-transformed (with SNR)
        fname = sprintf('sub-%03d_rhythm-%s_tone-%s_snr-%d-%d_FFT_dB',...
            subject, rhythm, tone, par.snr_bins_eeg(1), par.snr_bins_eeg(2)); 
        header_fft_db.name = fname; 
        CLW_save(fpath, header_fft_db, data_fft_db); 
    end
end

%% 

% take N best channels by considering general response magnitude as
% quality measure (this is not double dipping!)
n_best_chan = length(par.front_chan); 

% average across all conditions
amps_sum = mean(mean(amps_sum,3),2); 

[~,best_chan_idx] = sort(amps_sum,'descend'); 

chan_sel = {header_fft.chanlocs(best_chan_idx(1:n_best_chan)).labels}; 

% save data for selected channels
fpath = fullfile(par.deriv_path, ...
    sprintf('roi-best/sub-%03d',subject)); 

if ~isdir(fpath); mkdir(fpath); end 

% selected channel names 
if ~isempty(chan_sel)
    fname = sprintf('chan_names.csv'); 
    chan_table = cell2table(chan_sel', 'VariableNames',{'label'}); 
    writetable(chan_table,fullfile(fpath,fname)); 
end    

%% 

% update bad subjects table 
writetable(bad_subjects_table,fullfile(par.deriv_path,'bad_subjects.csv')); 









