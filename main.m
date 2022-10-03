clear all
addpath(genpath('./lib'))

par = getParams(); 

chan_sel_methods = {'front','best'}; 

importLetswave; 

%% 

% run cochlear model on the stimulus audio files
cochlear_model; 

% preprocess EEG data 
for iSub=1:length(par.subjects)
    preproc(par.subjects(iSub)); 
end

% get spectra and ERP for one rhythm cycle for ALL channels (this will be used 
% for topoplots and for selecting channels with greatest overall response)
for iSub=1:length(par.subjects)
    procAllChan(par.subjects(iSub)); 
end

% get spectra and ERP for one rhythm cycle for a selection of channels 
for iSub=1:length(par.subjects)
    for iSel=1:length(chan_sel_methods)
        procROI(par.subjects(iSub), chan_sel_methods{iSel})
    end
end

% extract magnitudes at frequencies of interest and save to text file that will
% be loaded in R
for iSel=1:length(chan_sel_methods)
    extractFeatures(chan_sel_methods{iSel})
end

% % plot and save a summary figure for each individual participant
% for iSub=1:length(par.subjects)
%     for iSel=1:length(chan_sel_methods)
%         plotSummary(par.subjects(iSub), chan_sel_methods{iSel})
%     end
% end

% plot and save a summary figure for grand average across all participants 
for iSel=1:length(chan_sel_methods)
    plotSummaryGrand(chan_sel_methods{iSel})
end

%%

% Next step is to run statistical analysis in R. Do this by manually executing 
% the `main.R` script in RStudio, which enables using a local project-specific 
% environment to ensure reproducibility. 

