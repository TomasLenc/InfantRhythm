function plotSummaryGrand(roi)

% prepare layout for fieldtrip topoplots 
importLetswave(-1); 
importFieldtrip; 

cfg  = [];
cfg.layout = 'GSN-HydroCel-129.sfp';
lay = ft_prepare_layout(cfg);

importFieldtrip(-1); 
importLetswave; 

par = getParams(); 

bad_subjects_table = readtable(fullfile(par.deriv_path,'bad_subjects.csv')); 
subjects2run = setdiff(par.subjects, bad_subjects_table.subject); 

%% prepare data


% sub, rhyth, tone
erp_all = cell(length(subjects2run), length(par.rhythms), length(par.tones)); 

% sub, rhyth, tone
data_fft_all = cell(length(subjects2run), length(par.rhythms), length(par.tones)); 

% sub, rhyth, tone
data_fft_noSNR_all = cell(length(subjects2run), length(par.rhythms), length(par.tones)); 

% sub, rhyth, tone, isMeterRel
amps_topo_all = cell(length(subjects2run), length(par.rhythms), length(par.tones), 2); 

% sub, rhyth, tone
zMeterRel_all = nan(length(subjects2run), length(par.rhythms), length(par.tones)); 


for iSub=1:length(subjects2run)

    subject = subjects2run(iSub); 
    fprintf('sub-%03d\n', subject); 
    
    for iRhythm=1:2

        for iTone=1:2

            tone = par.tones{iTone}; 
            rhythm = par.rhythms{iRhythm}; 

            % load time-domain data for selected channel 
            fpath = fullfile(par.deriv_path, ...
                sprintf('roi-%s/sub-%03d',roi,subject)); 
            fname = sprintf('sub-%03d_rhythm-%s_tone-%s_time',subject,rhythm,tone); 
            [header_time,data_time] = CLW_load(fullfile(fpath,fname)); 

            [t_erp, erp, sem_erp] = getCycleErp(data_time, 1/header_time.xstep, 2.4); 
            erp_all{iSub,iRhythm,iTone} = erp; 

            % FFT 
            fpath = fullfile(par.deriv_path, ...
                sprintf('roi-%s/sub-%03d',roi,subject)); 
            fname = sprintf('sub-%03d_rhythm-%s_tone-%s_snr-%d-%d_FFT',...
                subject,rhythm,tone,par.snr_bins_eeg(1), par.snr_bins_eeg(2)); 
            
            [header_fft,data_fft] = CLW_load(fullfile(fpath,fname)); 

            data_fft_all{iSub,iRhythm,iTone} = data_fft; 

            [zMeterRel,~,~,~,~,amps] = getZ(squeeze(data_fft), par.frex, ...
                par.idx_meterRel, par.idx_meterUnrel, 'xstep', header_fft.xstep); 
            zMeterRel_all(iSub,iRhythm,iTone) = zMeterRel; 
            
            % FFT no-SNR
            fpath = fullfile(par.deriv_path, ...
                sprintf('roi-%s/sub-%03d',roi,subject)); 
            fname = sprintf('sub-%03d_rhythm-%s_tone-%s_snr-0-0_FFT',...
                            subject,rhythm,tone); 
            [header_fft_noSNR,data_fft_noSNR] = CLW_load(fullfile(fpath,fname)); 
            data_fft_noSNR_all{iSub,iRhythm,iTone} = data_fft_noSNR; 
            
            % load FFT all channels (for topo) 
            fpath = fullfile(par.deriv_path, ...
                sprintf('roi-all/sub-%03d',subject)); 
            fname = sprintf('sub-%03d_rhythm-%s_tone-%s_snr-%d-%d_FFT',...
                subject,rhythm,tone,par.snr_bins_eeg(1), par.snr_bins_eeg(2)); 
            
            [header_fft_allChan,data_fft_allChan] = CLW_load(fullfile(fpath,fname)); 

            [~,~,~,ampsMeterRel,ampsMeterUnrel,amps] = getZ(squeeze(data_fft_allChan), par.frex, ...
                par.idx_meterRel, par.idx_meterUnrel, 'xstep', header_fft_allChan.xstep); 

            amps_topo_all{iSub,iRhythm,iTone,1} = ampsMeterRel'; 
            amps_topo_all{iSub,iRhythm,iTone,2} = ampsMeterUnrel'; 
            
            
        end
    end
end


% get limits for plotting
ampTopoLims = [0,0]; 
erpLim = 0; 
fftLims = [0,0]; 

for iRhythm=1:2
    for iTone=1:2
        
        erp = cat(1,erp_all{:,iRhythm,iTone}); 
        erp_mean = mean(erp, 1); 
        erp_sem = std(erp,[],1)/sqrt(size(erp,1)); 
        erpLim = max(erpLim, max(abs([erp_mean+erp_sem, erp_mean-erp_sem]))); 
        
        mX = cat(1,data_fft_all{:,iRhythm,iTone}); 
        [~,~,~,~,~,amps] = getZ(squeeze(mX), par.frex, ...
            par.idx_meterRel, par.idx_meterUnrel, 'xstep', header_fft.xstep); 
        amps_mean = mean(amps,1); 
        amps_sem = std(amps,[],1)/sqrt(size(amps,1)); 
        fftLims(2) = max(fftLims(2), max(abs(amps_mean+amps_sem))); 
        
        amps_topo_meterRel_grand = mean(cat(1,amps_topo_all{:,iRhythm,iTone,1}),1); 
        amps_topo_meterUnrel_grand = mean(cat(1,amps_topo_all{:,iRhythm,iTone,2}),1); 
        
        ampTopoLims(2) = max(ampTopoLims(2),...
            prctile([amps_topo_meterRel_grand,amps_topo_meterUnrel_grand],80)); 
        ampTopoLims(2) = ceil(ampTopoLims(2)*par.prec)/par.prec; 
        
    end
end

erpLim = ceil(erpLim*par.prec)/par.prec; 
fftLims(1) = floor(fftLims(1)*par.prec)/par.prec; 
fftLims(2) = ceil(fftLims(2)*par.prec)/par.prec; 
zLims = [-0.6, 1.4]; 

            
% channels used for topo highlighting
fpath = fullfile(par.deriv_path, ...
    sprintf('roi-%s/sub-%03d',roi,subject)); 
fname = sprintf('chan_names.csv'); 
if exist(fullfile(fpath,fname))
    chan_sel = readtable(fullfile(fpath,fname)); 
    chan_sel = chan_sel.label; 
else 
    chan_sel = []; 
end

%% 
    
% open figure and pack subplots
f = figure('color','white',...
           'name',sprintf('roi-%s',roi),...
           'position', [312 296 711 536] ); 
pnl = panel(f);

% rhythm
pnl.pack('v',2); 

% coch-time, mX, zscore
pnl(1).pack('h',[0.3, 0.4, 0.2]); 
pnl(2).pack('h',[0.3, 0.4, 0.2]); 

for iRhythm=1:2
    % time, mX
    for iPlotType=1:2
        % EEG low vs. EEG high 
        pnl(iRhythm,iPlotType).pack('v',2); 
    end
end

for iRhythm=1:2
    % EEG low vs. EEG high     
    for iTone=1:2
        pnl(iRhythm,1,iTone).pack({[0,0,1,1]}); 
        pnl(iRhythm,1,iTone).pack({[0,0.3,1,0.4]}); 
    end
end

for iRhythm=1:2
    % EEG low vs. EEG high     
    for iTone=1:2
        pnl(iRhythm,2,iTone).pack({[0,0,1,1]}); 
        % inset for ERP
        pnl(iRhythm,2,iTone).pack({[0.1, 0.9, 0.4, 0.3]}); 
        % inset for ampAround
        pnl(iRhythm,2,iTone).pack({[1.1, 0, 0.2, 0.3]}); 
        % inset for topo
        pnl(iRhythm,2,iTone).pack({[1.0, 0.6, 0.2, 0.6]}); 
        pnl(iRhythm,2,iTone).pack({[1.2, 0.6, 0.2, 0.6]}); 
    end
end

% zscore 
pnl(1,3).pack({[0.6, -0.1, 0.7, 0.25]}); 
pnl(2,3).pack({[0.6, 0.7, 0.7, 0.25]}); 

% pnl.select('all'); 


h = []; 

for iRhythm=1:2

    for iTone=1:2

            tone = par.tones{iTone}; 
            rhythm = par.rhythms{iRhythm}; 

            % coch
            % ----
            [s, fs_s, t_s] = loadStim(rhythm,tone,'cycle'); 
            s = s / max(s); 
            [coch, fs_coch, t_coch] = loadCochERP(rhythm,tone,'cycle');

            ax = pnl(iRhythm,1,iTone,2).select(); 
            hold(ax,'on'); 
            plot(ax, t_s, s, 'col',[.6,.6,.6]); 
            plot(ax, t_coch, coch, 'col',[0,0,0], 'linew',2); 
            ax.Visible = 'off'; 

            % ERP 
            % ---
            erp = cat(1,erp_all{:,iRhythm,iTone}); 
            erp_mean = mean(erp,1); 
            erp_sem = std(erp,[],1)/sqrt(size(erp,1)); 

            ax = pnl(iRhythm,2,iTone,2).select(); 
            plot(ax, t_s, s/max(s)*erpLim, 'col',[.8,.8,.8]); 
            hold(ax,'on')
            fill(ax, [t_erp,flip(t_erp)], ...
                 [erp_mean-erp_sem,flip(erp_mean+erp_sem)], ...
                 par.col_time_eeg(iTone,:), ...
                 'facealpha',0.3, ...
                 'LineStyle','none'); 
            h(iTone) = plot(ax, t_erp, erp_mean, ...
                'col',par.col_time_eeg(iTone,:), 'linew',1); 
            ax.XAxis.Visible = 'off'; 
            ax.YLim = [-erpLim,erpLim]; 
            ax.YTick = [-erpLim,erpLim]; 
            ax.YTickLabel = [-erpLim,erpLim]; 


            % spectrum
            % --------

            ax = pnl(iRhythm,2,iTone,1).select(); 

            mX = cat(1,data_fft_all{:,iRhythm,iTone}); 
            mX_mean = mean(mX, 1); 
            mX_sem = std(mX,[],1)/sqrt(size(mX,1)); 

            plotLwFFT(ax, header_fft, mX_mean, 'ci',mX_sem); 

            ax.YLim = [0,fftLims(2)]; 
            ax.YTick = [0,fftLims(2)]; 


            % amp-around
            % ----------
            mXnoSNR = squeeze(cat(1,data_fft_noSNR_all{:,iRhythm,iTone})); 

            [amp_around, amp_around_z] = getAmpAround( ...
                header_fft_noSNR, mXnoSNR, par.frex, par.snr_bins_eeg, par.amp_around_alpha); 

            ax = pnl(iRhythm,2,iTone,3).select(); 
            data2plot = squeeze(mean(amp_around,1)); 
            data2plot = data2plot-min(data2plot); 
            z2plot = squeeze(mean(amp_around_z,1)); 

            plotAmpAround(ax, data2plot, z2plot, par.snr_bins_eeg, par.prec); 


    end

    % zscores
    % -------
    z_meterRel = squeeze(zMeterRel_all(:,iRhythm,:)); 
    
    ax = pnl(iRhythm,3,1).select(); 
    hold on
   
    [mX_coch, xstep_coch] = loadCoch(rhythm,tone); 
    
    [z_meterRel_coch] = getZ(mX_coch, par.frex, ...
        par.idx_meterRel, par.idx_meterUnrel, 'xstep', xstep_coch); 

    % test each task against cochlear model 
    plotCondDiffPaired(ax, z_meterRel, ...
        par.col_neutral, par.col_z_eeg, 'mu', z_meterRel_coch);
    
    ax.YLim = zLims; 

    
end


%% margins and labels (main figure) 

% [left bottom right top]
pnl.de.margin = [10,4,0,10]; 

pnl(1,3).marginleft = 25; 
pnl(2,3).marginleft = 25; 

% yticklabels are overlapping (2.4 0), fix it
pnl(2).margintop = 15; 

pnl.margin = [10,5,6,10];  

% add legend 
pnl(iRhythm,3,1).select(); 

l = legend(h, {par.tones{1}',par.tones{2}}); 
l.Box = 'off'; 
l.FontSize = par.fontsize; 
l.Position = [0.8486 0.1999 0.1058 0.0822]; 

pnl.fontsize = par.fontsize; 


%% save (main figure)

fpath = fullfile(par.figures_path, ...
                sprintf('roi-%s',roi)); 

if ~isfolder(fpath)
    mkdir(fpath); 
end

fname = sprintf('roi-%s_summaryGrand',roi); 

saveas(f, fullfile(fpath,[fname,'.svg'])); 
print(f, '-dpng', '-painters', '-r300', fullfile(fpath,fname)); 

close(f); 



     
%% figure topo 


% Topoplots never render well in a composite figure, and so it's better to just
% plot the separately, save as raster and then paste manually in inkscape. 

% open figure and pack subplots
f_topo = figure('color','white',...
           'name',sprintf('roi-%s',roi),...
           'position', [810 565 708 315] ); 
pnl_topo = panel(f_topo);

% rhythm
pnl_topo.pack('v',2); 

% tone
pnl_topo(1).pack('h',2); 
pnl_topo(2).pack('h',2); 

for iRhythm=1:2
    for iTone=1:2
        pnl_topo(iRhythm,iTone).pack('h',2); 
    end
end

% pnl_topo.select('all')


for iRhythm=1:2

    for iTone=1:2

        tone = par.tones{iTone}; 
        rhythm = par.rhythms{iRhythm}; 

        % topo
        % ----
        amps_topo_meterRel_mean = mean(cat(1,amps_topo_all{:,iRhythm,iTone,1}),1); 
        amps_topo_meterUnrel_mean = mean(cat(1,amps_topo_all{:,iRhythm,iTone,2}),1); 
         
        importLetswave(-1); 
        importFieldtrip; 
        
        [cfg,d] = prepareTopoCfg(amps_topo_meterRel_mean, ...
                                {header_fft_allChan.chanlocs.labels}, ...
                                ampTopoLims, lay); 
        ax = axes(f_topo); 
        ft_topoplotER(cfg,d);
        colormap(parula) 
        pnl_topo(iRhythm,iTone,1).select(ax); 
        
        [cfg,d] = prepareTopoCfg(amps_topo_meterUnrel_mean, ...
                                {header_fft_allChan.chanlocs.labels}, ...
                                ampTopoLims, lay); 
        ax = axes(f_topo); 
        ft_topoplotER(cfg,d);
        colormap(parula) 
        pnl_topo(iRhythm,iTone,2).select(ax); 

        importFieldtrip(-1); 
        importLetswave; 

    end
end


pnl_topo.de.margin = 0; 
pnl_topo.margin = [0,0,0,10]; 

pnl_topo(1).title(par.rhythms(1)); 
pnl_topo(2).title(par.rhythms(2)); 
pnl_topo(1,1).title(par.tones(1)); 
pnl_topo(1,2).title(par.tones(2)); 
pnl_topo(1,1,1).title('meter-rel'); 
pnl_topo(1,1,2).title('meter-unrel'); 


fname = sprintf('roi-%s_summaryGrand-topo',roi); 
print(f_topo, '-dpng', '-painters', '-r600', fullfile(fpath,fname)); 
close(f_topo); 


f_cbar = figure('color','white', 'pos',[704 542 536 112]); 
axes();
cbar = colorbar();
cbar.Ticks = [0,1]; 
cbar.TickLabels = ampTopoLims; 
cbar.Label.String = 'amplitude'; 
fname = sprintf('roi-%s_summaryGrand-topo-cbar',roi); 
print(f_cbar, '-dpng', '-painters', '-r600', fullfile(fpath,fname)); 
close(f_cbar); 



















