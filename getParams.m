function params=getParams()

% ========== manually update these paths before running the analysis ==========
% path the the experiment root folder 
experiment_path = '/datadisk/Dropbox/projects/XPInfant'; 

% path the letswave6
letswave_path = '/datadisk/Dropbox/libraries/matlab/_neuro/letswave6-master';

% path to fieldtrip
fieldtrip_path = '/home/tomo/Documents/MATLAB/fieldtrip';
% ============================================================================

source_path = fullfile(experiment_path,'source'); 
stim_path = fullfile(experiment_path,'stimuli'); 

deriv_path = fullfile(experiment_path,'derivatives'); 
if ~isdir(deriv_path)
    mkdir(deriv_path)
end
figures_path = fullfile(experiment_path,'figures'); 
if ~isdir(figures_path)
    mkdir(figures_path)
end
preproc_path = fullfile(experiment_path,'derivatives/preprocessed_AB'); 
if ~isdir(preproc_path)
    mkdir(preproc_path)
end
features_path = fullfile(experiment_path,'derivatives/features'); 
if ~isdir(features_path)
    mkdir(features_path)
end

%%

subjects = [1:20]; 

tones = {'low','high'}; 
rhythms = {'unsync','sync'}; 


trial_duration = 60;

cycle_duration = 2.4; 


% mastoids 
ref_chans = {'E57','E100'}; 

% set this to true if you want the preprocessed trials to be chunked before
% averaging (this might help to increase the SNR?)
preproc_chunk_do = true; 
preproc_chunk_onset = 2.4; 
preproc_chunk_dur = 26.4; 


rois = {'best','front'}; 

n_chan = 92; 

front_chan = {'E40','E41','E34','E35','E36','E27','E28','E29','E24','E20','E19','E46',...
        'E4','E118','E123','E124','E111','E116','E117','E109','E110','E111','E102','E103','E104'}; 


% number of bootstrap samples to estimate confidence intervals 
nBoot = 200; 

%% frequencies of interest

% maximum limit for plotting frequency axis
maxfreqlim = 5.5; % 5.5

% original Sylvie's selection: 
frex = 1/2.4 * [2:12]; 
idx_meterRel = [2,5,8,11]; 

% % without 0.416 and 0.833 (lowest two frex): 
% frex = 1/2.4 * [3:12]; 
% idx_meterRel = [1,4,7,10]; 

% % without 0.416, 0.833 and 5 Hz: 
% frex = 1/2.4 * [3:11]; 
% idx_meterRel = [1,4,7]; 


idx_meterUnrel = setdiff([1:length(frex)],idx_meterRel); 

% this is for extraction to csv file (we'll process this further in R)
frex_to_extract = 1/2.4 * [1:72]; % take all frequencies up to 30Hz


%% FFT 

snr_bins_eeg = [2,5]; 

snr_bins_coch = [2,5]; 

amp_around_alpha = 0.001; 




%% plotting 

col_meterRel = [222 45 38]/255; 
col_meterUnrel = [49, 130, 189]/255; 
col_neutral = repmat(0.5,1,3); 

col_time_eeg = [138, 54, 129; 
                16, 125, 52]/255; 
col_z_eeg = [brighten(col_time_eeg(1,:), 0.6), 
             brighten(col_time_eeg(2,:), 0.6)]; 

col_time_sound = [0,0,0]; 
col_z_sound = [0,0,0]; 

linew_fft = 1.7; 

fontsize = 12; 


prec = 100; 


%% return structure 

params = []; 

w = whos;
for a = 1:length(w) 
    params.(w(a).name) = eval(w(a).name); 
end


