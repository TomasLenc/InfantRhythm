function plotSummary(subject, chan_sel_method)

% clear all
% close all


% prepare layout for fieldtrip topoplots 
importLetswave(-1); 
importFieldtrip; 

cfg  = [];
cfg.layout = 'GSN-HydroCel-129.sfp';
lay = ft_prepare_layout(cfg);

importFieldtrip(-1); 
importLetswave; 

%% 



% iSub=1; 
% iTone = 1; 
% iRhythm = 1; 
% 
% subject = subjects{iSub}; 
% tone = tones{iTone}; 
% rhythm = rhythms{iRhythm}; 

par = getParams(); 


%% 
    


% open figure and pack subplots
f = figure('color','white',...
           'name',sprintf('sub-%03d roi-%s',subject,chan_sel_method),...
           'position', [743 80 1089 761]); 
pnl = panel(f);

% rhythm
pnl.pack('v',2); 

% time, mX, zscore
pnl(1).pack('h',[0.4,0.3,0.2]); 
pnl(2).pack('h',[0.4,0.3,0.2]); 

for iRhythm=1:2
    % time, mX
    for iPlotType=1:2
        % coch vs. EEG low vs. EEG high 
        pnl(iRhythm,iPlotType).pack('v',3); 
    end
end

for iRhythm=1:2
    % EEG low vs. EEG high     
    for iTone=1:2
        pnl(iRhythm,2,iTone+1).pack({[0,0,1,1]}); 
        % inset for ampAround
        pnl(iRhythm,2,iTone+1).pack({[1.1, 0, 0.2, 0.3]}); 
        % inset for topo
        pnl(iRhythm,2,iTone+1).pack({[1.0, 0.6, 0.4, 0.6]}); 
    end
end

% zscore 
pnl(1,3).pack({[0.6, -0.1, 0.7, 0.25]}); 
pnl(2,3).pack({[0.6, 0.7, 0.7, 0.25]}); 

%     pnl.select('all'); 
%     

%% prepare limits 

ampSumTopoLims = [0,0]; 
erpLim = 0; 
fftLims = [0,0]; 

for iRhythm=1:2

    for iTone=1:2

        tone = par.tones{iTone}; 
        rhythm = par.rhythms{iRhythm}; 

        % load FFT all channels (for topo) 
        fpath = fullfile(par.deriv_path, ...
            sprintf('roi-all/sub-%03d',subject)); 
        fname = sprintf('sub-%03d_rhythm-%s_tone-%s_snr-%d-%d_FFT',...
            subject,rhythm,tone,par.snr_bins_eeg(1), par.snr_bins_eeg(2)); 
        [header_fft_all,data_fft_all] = CLW_load(fullfile(fpath,fname)); 

        [~,~,~,amps,~] = getZ(squeeze(data_fft_all), par.frex, ...
            par.idx_meterRel, par.idx_meterUnrel, 'xstep', header_fft_all.xstep); 
        amp_sum = sum(amps,2); 
        ampSumTopoLims(2) = max(ampSumTopoLims(2), prctile(amp_sum,80)); 

        % laod time-domain data for selected channel 
        fpath = fullfile(par.deriv_path, ...
            sprintf('roi-%s/sub-%03d',chan_sel_method,subject)); 

        fname = sprintf('sub-%03d_rhythm-%s_tone-%s_time',subject,rhythm,tone); 
        [header_time,data_time] = CLW_load(fullfile(fpath,fname)); 

        [t_erp, erp, sem_erp] = getCycleErp(data_time, 1/header_time.xstep, 2.4); 
        erpLim = max(erpLim, max(abs(erp))); 

        % laod freq-domain data for selected channel 
        fname = sprintf('sub-%03d_rhythm-%s_tone-%s_snr-%d-%d_FFT',...
            subject,rhythm,tone,par.snr_bins_eeg(1), par.snr_bins_eeg(2)); 
        [header_fft,data_fft] = CLW_load(fullfile(fpath,fname)); 

        [~,~,~,amps,~] = getZ(squeeze(data_fft), par.frex, ...
            par.idx_meterRel, par.idx_meterUnrel, 'xstep', header_fft.xstep); 
        fftLims(1) = min(fftLims(1), min(squeeze(data_fft))); 
        fftLims(2) = max(fftLims(2), max(abs(amps))); 

    end
end

ampSumTopoLims(2) = ceil(ampSumTopoLims(2)*par.prec)/par.prec; 
erpLim = ceil(erpLim*par.prec)/par.prec; 
fftLims(1) = floor(fftLims(1)*par.prec)/par.prec; 
fftLims(2) = ceil(fftLims(2)*par.prec)/par.prec; 


%%
h = []; 

for iRhythm=1:2

    for iTone=1:2

        tone = par.tones{iTone}; 
        rhythm = par.rhythms{iRhythm}; 

        % load data
        % ---------

        % load FFT all channels (for topo) 
        fpath = fullfile(par.deriv_path, ...
            sprintf('roi-all/sub-%03d',subject)); 
        fname = sprintf('sub-%03d_rhythm-%s_tone-%s_snr-%d-%d_FFT',...
            subject,rhythm,tone,par.snr_bins_eeg(1), par.snr_bins_eeg(2)); 
        [header_fft_all,data_fft_all] = CLW_load(fullfile(fpath,fname)); 

        for iChan=1:length(header_fft_all.chanlocs)
            idx = find(strcmp(lay.label,header_fft_all.chanlocs(iChan).labels)); 
            header_fft_all.chanlocs(iChan).x = lay.pos(idx,1); 
            header_fft_all.chanlocs(iChan).y = lay.pos(idx,2); 
        end

        % laod data for selected channel 
        fpath = fullfile(par.deriv_path, ...
            sprintf('roi-%s/sub-%03d',chan_sel_method,subject)); 

        % channels used 
        fname = sprintf('chan_names.csv'); 
        if exist(fullfile(fpath,fname))
            chan_sel = readtable(fullfile(fpath,fname)); 
            chan_sel = chan_sel.label; 
        else 
            chan_sel = []; 
        end

        % time-domain 
        fname = sprintf('sub-%03d_rhythm-%s_tone-%s_time',subject,rhythm,tone); 
        [header_time,data_time] = CLW_load(fullfile(fpath,fname)); 

        % FFT 
        fname = sprintf('sub-%03d_rhythm-%s_tone-%s_snr-%d-%d_FFT',...
            subject,rhythm,tone,par.snr_bins_eeg(1), par.snr_bins_eeg(2)); 
        [header_fft,data_fft] = CLW_load(fullfile(fpath,fname)); 

        % FFT no-SNR
        fname = sprintf('sub-%03d_rhythm-%s_tone-%s_snr-0-0_FFT',...
                        subject,rhythm,tone); 
        [header_fft_noSNR,data_fft_noSNR] = CLW_load(fullfile(fpath,fname)); 


        % ERP 
        % ---

        [t_erp, erp, sem_erp] = getCycleErp(data_time, 1/header_time.xstep, 2.4); 

        [s, fs_s, t_s] = loadStim(rhythm,tone,'cycle'); 
        s = s/max(s) * erpLim; 

        ax = pnl(iRhythm,1,iTone+1).select(); 
        plot(ax, t_s, s, 'col',[.8,.8,.8]); 
        hold(ax,'on')
        fill(ax, [t_erp,flip(t_erp)], ...
             [erp-sem_erp,flip(erp+sem_erp)], ...
             par.col_time_eeg(iTone,:), ...
             'facealpha',0.3, ...
             'LineStyle','none'); 
        h(iTone) = plot(ax, t_erp, erp, 'col',par.col_time_eeg(iTone,:), 'linew',2); 
        ax.XAxis.Visible = 'off'; 
        ax.YLim = [-erpLim,erpLim]; 

        
        % spectrum
        % --------

        ax = pnl(iRhythm,2,iTone+1,1).select(); 

        plotLwFFT(ax, header_fft, data_fft); 

        ax.YLim = [0,fftLims(2)]; 
        ax.YTick = [0,fftLims(2)]; 

        
        % amp-around
        % ----------

        mXnoSNR = squeeze(data_fft_noSNR)'; 

        [amp_around, amp_around_z] = getAmpAround( ...
            header_fft_noSNR, mXnoSNR, par.frex, par.snr_bins_eeg, par.amp_around_alpha); 

        ax = pnl(iRhythm,2,1+iTone,2).select(); 
        data2plot = squeeze(mean(amp_around,1)); 
        z2plot = squeeze(mean(amp_around_z,1)); 
        
        plotAmpAround(ax, data2plot, z2plot, par.snr_bins_eeg, par.prec); 

        
        % topo
        % ----

        importLetswave(-1); 
        importFieldtrip; 

        [~,~,~,~,~,amps] = getZ(squeeze(data_fft_all), par.frex, ...
            par.idx_meterRel, par.idx_meterUnrel, 'xstep', header_fft_all.xstep); 

        amps_sum = sum(amps,2);

        d = struct;
        d.avg = amps_sum; 
        d.time = 1;
        d.label = {header_fft_all.chanlocs.labels};
        d.dimord = 'chan_time';

        cfg = [];                            
        cfg.xlim = [1 1];   % time to plot             
        cfg.zlim = ampSumTopoLims;  % amplitude range              
        cfg.layout = lay;     
        cfg.marker = 'off'; 
        cfg.comment = ' ';   
        cfg.commentpos = 'title'; 
        cfg.gridscale = 500;   
        cfg.style = 'straight';
        cfg.highlight = 'on'; 
        cfg.highlightsymbol = '.'; 
        cfg.highlightsize = 5; 
        cfg.highlightcolor = 'k'; 
        cfg.highlightchannel = chan_sel; 

        ax = axes; 
        ft_topoplotER(cfg,d);
        colormap(parula) 

        pnl(iRhythm,2,1+iTone,3).select(ax); 

        importFieldtrip(-1); 
        importLetswave; 


        % zscores
        % -------

        % calculate amplitudes and zscores 
        [z_meterRel,z_meterUnrel] = getZ(squeeze(data_fft), par.frex, ...
            par.idx_meterRel, par.idx_meterUnrel, 'xstep',header_fft.xstep); 

        % plot 
        pnl(iRhythm,3,1).select(); 

        hold on
        plot([0,1],[z_meterRel,z_meterUnrel],'-o',...
                    'Color',par.col_z_eeg(iTone,:), ...
                    'MarkerEdgeColor',par.col_z_eeg(iTone,:),...
                    'MarkerSize',8,...
                    'LineWidth',2)
        xlim([-0.5,1.5])
        ylim([-1.2,1.2])
        ax = gca; 
        ax.Color = 'none'; 
        ax.XTick = []; 
        ax.XAxisLocation = 'origin'; 
        ax.XTickLabel = []; 
        ax.YTick = [-1,1]; 
        ax.FontSize = par.fontsize; 
        ax.LineWidth = 2;


    end
end



% ------------
% MARGINS and LABELS
% ------------

% [left bottom right top]
pnl.de.margin = [10,0,0,6]; 

pnl(1,3).marginleft = 15; 

% yticklabels are overlapping (2.4 0), fix it
pnl(2).margintop = 20; 

pnl.margin = [17,15,15,8];  

% add legend 
pnl(iRhythm,3,1).select(); 

l = legend(h, {par.tones{1}',par.tones{2}}); 
l.Box = 'off'; 
l.FontSize = par.fontsize; 
l.Position = [0.8486 0.1999 0.1058 0.0822]; 

pnl.fontsize = par.fontsize; 

% % add title 
% pnl(1,2,1).select(); 
% tit = title(sprintf('sub-%03d chan-%s',subject,chan_sel_method), 'Interpreter','none'); 
% tit.FontSize = par.fontsize; 
% tit.FontWeight = 'bold'; 



% ------------
% SAVE
% ------------

fpath = fullfile(par.figures_path, ...
                sprintf('roi-%s',chan_sel_method)); 

if ~isfolder(fpath)
    mkdir(fpath); 
end

fname = sprintf('sub-%03d_roi-%s',subject,chan_sel_method); 

saveas(f, fullfile(fpath,[fname,'.fig'])); 
print(f, '-dpng', '-painters', fullfile(fpath,fname)); 

close(f); 



    
    

















