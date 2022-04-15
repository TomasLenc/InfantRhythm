function [amp_around,amp_around_z,amp_around_signif]=getAmpAround( ...
    header, mXnoSNR, frex, snr_bins, alpha)

% number of frequency bins on each side of the target frequency bin
n_sidebins = snr_bins(2); 
sidebin_idx = [-n_sidebins:+n_sidebins]; 

% so we need to allocate 
n_bins = 2*n_sidebins+1;

%% 

frex_idx = round(frex/header.xstep)+1; 

% allocate [channel,bin]
amp_around = zeros(size(mXnoSNR,1), n_bins); 

% go across harmonics, take a segment of the spectra around, and sum across
% harmonics (will use this for z-score SNR) 
for fi=1:length(frex_idx)
    amp_around = amp_around + mXnoSNR(:,frex_idx(fi)+sidebin_idx); 
end

% calculate significance for each contact 
signal = amp_around(:,n_sidebins+1); 

noise_bins = [-snr_bins(2):-snr_bins(1), ...
              snr_bins(1):snr_bins(2)]; 

noise = amp_around(:,[n_sidebins+1 + noise_bins]);  
noise_mean = mean(noise,2); 
noise_sd = std(noise,[],2); 

amp_around_z = (signal-noise_mean)./noise_sd; 

p = 1-normcdf(amp_around_z); 
amp_around_signif = p < alpha; 

