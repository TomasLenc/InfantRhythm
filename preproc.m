function preproc(subject)

par = getParams(); 

rhythms = {'unsync','sync'}; 
tones = {'low','high'}; 

for iRhythm=1:length(rhythms)
    
    for iTone=1:length(tones)
        
        d = dir(fullfile(par.source_path, ...
            sprintf('%s %s P%03d*.mat', tones{iTone}, rhythms{iRhythm}, subject))); 
        
        if length(d) > 1
            error('too many source files matching sub-%03d %s %s', ...
                subject, tones{iTone}, rhythms{iRhythm}); 
        else
            fprintf('\nsub-%03d %s %s\n', subject, tones{iTone}, rhythms{iRhythm}); 
        end
        % check if Kingswood 
        if strfind(d.name, 'K.mat')
            isKingswood = true; 
        else
            isKingswood = false; 
        end

        % load data
        [header, data] = CLW_load(fullfile(d.folder, d.name));     
        
        % segment -1 to 60s
        trig_code = unique({header.events.code}); 
        assert(length(trig_code)==1)
        [header, data] = RLW_segmentation(header, data, trig_code, ...
            'x_start', -0.5, 'x_duration', 61); 
    
        % filter
        [header, data] = RLW_butterworth_filter(header, data, ...
            'filter_type', 'bandpass', ...
            'low_cutoff', 0.5, ...
            'high_cutoff', 30, ...
            'filter_order', 4); 
        
        % decimate factor 4
        [header, data] = RLW_downsample(header, data, ...
            'x_downsample_ratio', 4); 
        
        % remove DC
        [header, data] = RLW_dc_removal(header, data); 
        
        % segemnt from 0 to 60s
        [header, data] = RLW_segmentation(header, data, trig_code, ...
            'x_start', 0, 'x_duration', 60); 
            
        % remove 3 outer rings of channels 
        chan2rm = table2cell(readtable(fullfile(par.experiment_path,'code','chan2rm'), ...
            'ReadVariableNames',false)); 
        chan2rm = cellfun(@(x) ['E',num2str(x)], chan2rm, 'uni',0); 
        chan2keep = setdiff({header.chanlocs.labels}, chan2rm); 
        chan2keep = natsort(chan2keep); 
        [header, data] = RLW_arrange_channels(header, data, chan2keep); 
        
        % artifact blocking 
        data_tmp = squeeze(permute(data,[6,5,4,3,1,2])); 
        data_concat = reshape(data_tmp, [], length(chan2keep))'; 

        cfg = [];
        cfg.Approach = 'Window'; 
        cfg.Threshold = 100; % the threshold for detecting artifact
        cfg.Fs = 1/header.xstep; % sampling frequency
        cfg.WindowSize = 10; % unit in second
        cfg.InData = data_concat;  % the EEG data
        
        cfg = Run_AB(cfg);

        n_epochs = size(data,1); 
        ab_data = reshape(cfg.OutData', header.datasize(end), n_epochs, length(chan2keep));
        ab_data = permute(ab_data, [2,3,1]); 
        out_data = []; 
        out_data(:,:,1,1,1,:) = ab_data;
        
        % save data
        out_path = par.preproc_path; 
        if ~isdir(out_path)
            mkdir(out_path)
        end
        out_header = header; 
        out_header.name = ['AB ', header.name]; 
        CLW_save(out_path, out_header, out_data); 
        
    end
end
