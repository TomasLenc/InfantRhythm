function plotSummaryGrand(roi, varargin)

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

fft_type_eeg = 'amp';
fft_type_coch = 'amp';
if any(strcmpi(varargin, 'fft_type_eeg'))
    fft_type_eeg = varargin{find(strcmpi(varargin, 'fft_type_eeg')) + 1}; 
end
if any(strcmpi(varargin, 'fft_type_coch'))
    fft_type_coch = varargin{find(strcmpi(varargin, 'fft_type_coch')) + 1}; 
end

%% prepare data


% sub, rhyth, tone
erp_all = cell(length(subjects2run), length(par.rhythms), length(par.tones)); 

% sub, rhyth, tone
data_fft_all = cell(length(subjects2run), length(par.rhythms), length(par.tones)); 

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

            % load FFT averaged across roi channels
            fpath = fullfile(par.deriv_path, ...
                sprintf('roi-%s/sub-%03d',roi,subject)); 
            if strcmpi(fft_type_eeg, 'db')
                suffix = '_dB'; 
            else
                suffix = '';
            end
            fname = sprintf(...
                'sub-%03d_rhythm-%s_tone-%s_snr-%d-%d_FFT%s',...
                subject, rhythm, tone, ...
                par.snr_bins_eeg(1), par.snr_bins_eeg(2), ...
                suffix...
                ); 
            
            [header_fft, data_fft] = CLW_load(fullfile(fpath,fname)); 

            data_fft_all{iSub,iRhythm,iTone} = data_fft; 

            [zMeterRel,~,~,~,~,amps] = getZ(...
                squeeze(data_fft), par.frex, ...
                par.idx_meterRel, par.idx_meterUnrel, ...
                'xstep', header_fft.xstep...
                ); 
            zMeterRel_all(iSub,iRhythm,iTone) = zMeterRel; 
                        
            % load FFT for each channel (for topo) 
            fpath = fullfile(par.deriv_path, ...
                sprintf('roi-all/sub-%03d',subject)); 
            if strcmpi(fft_type_eeg, 'db')
                suffix = '_dB'; 
            else
                suffix = '';
            end
            fname = sprintf(...
                'sub-%03d_rhythm-%s_tone-%s_snr-%d-%d_FFT%s',...
                subject, rhythm, tone, ...
                par.snr_bins_eeg(1), par.snr_bins_eeg(2), ...
                suffix...
                ); 
            
            [header_fft_allChan, data_fft_allChan] = CLW_load(fullfile(fpath,fname)); 

            [~,~,~, ampsMeterRel, ampsMeterUnrel, amps] = getZ(...
                squeeze(data_fft_allChan), par.frex, ...
                par.idx_meterRel, par.idx_meterUnrel, ...
                'xstep', header_fft_allChan.xstep...
                ); 

            amps_topo_all{iSub,iRhythm,iTone,1} = ampsMeterRel'; 
            amps_topo_all{iSub,iRhythm,iTone,2} = ampsMeterUnrel'; 
            
            
        end
    end
end


% get limits for plotting
ampTopoLims = [0,0]; 
erpLim = 0; 
fftLims = {[0,0], [0,0]}; 

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
        fftLims{iRhythm}(2) = max(fftLims{iRhythm}(2), max(abs(amps_mean+amps_sem))); 
        
        amps_topo_meterRel_grand = mean(cat(1,amps_topo_all{:,iRhythm,iTone,1}),1); 
        amps_topo_meterUnrel_grand = mean(cat(1,amps_topo_all{:,iRhythm,iTone,2}),1); 
        
        ampTopoLims(2) = max(ampTopoLims(2),...
            prctile([amps_topo_meterRel_grand,amps_topo_meterUnrel_grand],80)); 
        ampTopoLims(2) = ceil(ampTopoLims(2)*par.prec)/par.prec; 
        
    end
    fftLims{iRhythm}(1) = floor(fftLims{iRhythm}(1)*par.prec)/par.prec; 
    fftLims{iRhythm}(2) = ceil(fftLims{iRhythm}(2)*par.prec)/par.prec; 
end

erpLim = ceil(erpLim*par.prec)/par.prec; 
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
           'position',[302 585 703 223]); 
pnl = panel(f);

% rhythm
pnl.pack('v',2); 

% low high 
pnl(1).pack('h',2); 
pnl(2).pack('h',2); 

for iRhythm=1:2
    for iTone=1:2
        % coch, EEG
        pnl(iRhythm,iTone).pack('h',2); 
    end
end

for iRhythm=1:2
    for iTone=1:2
        for iType=1:2
            pnl(iRhythm,iTone,iType).pack({[0,0,1,1]}); 
            % inset for ERP
            pnl(iRhythm,iTone,iType).pack({[0.4, 1.05, 0.6, 0.3]}); 
        end
    end
end

% pnl.select('all'); 


h = []; 

for iRhythm=1:2

    for iTone=1:2

            tone = par.tones{iTone}; 
            rhythm = par.rhythms{iRhythm}; 

            % coch
            % ----
            [s, ~, t_s] = loadStim(rhythm, tone, 'cycle'); 
            s = s / max(s); 
            [coch, ~, t_coch] = loadCochERP(rhythm, tone, 'cycle');

            ax = pnl(iRhythm,iTone,1,2).select(); 
            hold(ax,'on'); 
            plot(ax, t_s, s, 'col',[.8,.8,.8]); 
            plot(ax, t_coch, coch, 'col',[0,0,0], 'linew',1); 
            ax.Visible = 'off'; 
            
            [header_coch, data_coch] = loadCochFFT(...
                rhythm, tone, fft_type_coch, ...
                'maxfreqlim', par.maxfreqlim...
                );
             
            ax = pnl(iRhythm,iTone,1,1).select(); 
            plotLwFFT(ax, header_coch, data_coch); 
            ax.YLim = [0,Inf]; 
            ax.YTick = []; 
            
            % ERP 
            % ---
            erp = cat(1,erp_all{:,iRhythm,iTone}); 
            erp_mean = mean(erp,1); 
            erp_sem = std(erp,[],1)/sqrt(size(erp,1)); 

            ax = pnl(iRhythm,iTone,2,2).select(); 
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

            ax = pnl(iRhythm,iTone,2,1).select(); 

            mX = cat(1,data_fft_all{:,iRhythm,iTone}); 
            mX_mean = mean(mX, 1); 
            mX_sem = std(mX,[],1)/sqrt(size(mX,1)); 

%             plotLwFFT(ax, header_fft, mX_mean, 'ci',mX_sem); 
            plotLwFFT(ax, header_fft, mX_mean); 

            ax.YLim = [0,fftLims{iRhythm}(2)]; 
            ax.YTick = [0,fftLims{iRhythm}(2)]; 

    end
    
end
    

% margins and labels
% ------------------

% [left bottom right top]
pnl.de.margin = [10,4,0,10]; 

pnl.margin = [5,5,5,10];  

pnl.fontsize = 10; 


%% save (main figure)

fpath = fullfile(par.figures_path, ...
                sprintf('roi-%s',roi)); 

if ~isfolder(fpath)
    mkdir(fpath); 
end

fname = sprintf('roi-%s_summaryGrand_fftTypeEEG-%s_fftTypeCoch-%s', ...
                roi, fft_type_eeg, fft_type_coch); 

saveas(f, fullfile(fpath,[fname,'.svg'])); 
print(f, '-dpng', '-painters', '-r600', fullfile(fpath,fname)); 

close(f); 



     
%% figure topo 

% Topoplots never render well in a composite figure, and so it's better to just
% plot the separately, save as raster and then paste manually in
% inkscape...sorry

if ~any(strcmpi(varargin, 'skip_topo'))

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
        f_topo = figure('color','none',...
                   'name',sprintf('roi-%s',roi),...
                   'position', [1357 21 524 241] ); 
        pnl_topo = panel(f_topo);                            
        ax = axes(f_topo); 
        ft_topoplotER(cfg,d);
        colormap(parula) 
        pnl_topo().select(ax); 
        pnl_topo.de.margin = 0; 
        pnl_topo.margin = [1,1,1,1]; 
        fname = sprintf(...
            'roi-%s_rhythm-%s_tone-%s_fftTypeEEG-%s_meter-rel_summaryGrand-topo', ...
            roi, rhythm, tone, fft_type_eeg...
            ); 
        export_fig(f_topo, fullfile(fpath,fname), '-dpng', '-transparent', '-r300');
        close(f_topo); 
        
        [cfg,d] = prepareTopoCfg(amps_topo_meterUnrel_mean, ...
                                {header_fft_allChan.chanlocs.labels}, ...
                                ampTopoLims, lay); 
        f_topo = figure('color','none',...
                   'name',sprintf('roi-%s',roi),...
                   'position',[1357 21 524 241] ); 
        pnl_topo = panel(f_topo);                            
        ax = axes(f_topo); 
        ft_topoplotER(cfg,d);
        colormap(parula) 
        pnl_topo().select(ax); 
        pnl_topo.de.margin = 0; 
        pnl_topo.margin = [1,1,1,1]; 
        fname = sprintf(...
            'roi-%s_rhythm-%s_tone-%s_fftTypeEEG-%s_meter-unrel_summaryGrand-topo', ...
            roi, rhythm, tone, fft_type_eeg...
            ); 
        export_fig(f_topo, fullfile(fpath,fname), '-dpng', '-transparent', '-r300');
        close(f_topo); 
        
        importFieldtrip(-1); 
        importLetswave; 
    end
end

f_cbar = figure('color','white', 'pos',[704 542 536 112]); 
axes();
cbar = colorbar();
cbar.Ticks = [0,1]; 
cbar.TickLabels = ampTopoLims; 
cbar.Label.String = fft_type_eeg; 
fname = sprintf('roi-%s_fftTypeEEG-%s_summaryGrand-topo-cbar', ...
    roi, fft_type_eeg); 
print(f_cbar, '-dpng', '-painters', '-r600', fullfile(fpath,fname)); 
close(f_cbar); 


end















