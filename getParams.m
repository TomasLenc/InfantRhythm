function params=getParams()

[~, hostname] = system('hostname');
hostname = deblank(hostname); 

if strcmp(hostname, 'tux')
    deriv_path = '/datadisk/Dropbox/projects/XPInfant/derivatives'; 
    stim_path = '/datadisk/Dropbox/projects/XPInfant/stimuli'; 
    figures_path = '/datadisk/Dropbox/projects/XPInfant/figures'; 
    preproc_path = '/datadisk/Dropbox/projects/XPInfant/derivatives/preprocessed_AB'; 
    features_path = '/datadisk/Dropbox/projects/XPInfant/derivatives/features'; 
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

col_time_eeg = [64,0,166; 
                201, 96, 26]/255; 
col_z_eeg = [brighten(col_time_eeg(1,:), 0.6), 
             brighten(col_time_eeg(2,:), 0.6)]; 

col_time_sound = [0,0,0]; 
col_z_sound = [0,0,0]; 

linew_fft = 1.7; 

fontsize = 12; 


prec = 100; 


%% return structure 

% I'm lazy so here I take the workspace of this function and shove all the
% variables into a structure.  

% A little bonus is that the variables are alphabetically ordered :) 

params = []; 

w = whos;
for a = 1:length(w) 
    params.(w(a).name) = eval(w(a).name); 
end


