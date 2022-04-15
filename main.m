clear all
addpath(genpath('./lib'))

par = getParams(); 

chan_sel_methods = {'best','front'}; 

importLetswave; 

for iSub=1:length(par.subjects)
    procAllChan(par.subjects(iSub)); 
end

for iSub=1:length(par.subjects)
    for iSel=1:length(chan_sel_methods)
        procROI(par.subjects(iSub), chan_sel_methods{iSel})
    end
end

for iSel=1:length(chan_sel_methods)
    extractFeatures(chan_sel_methods{iSel})
end

for iSub=1:length(par.subjects)
    for iSel=1:length(chan_sel_methods)
        plotSummary(par.subjects(iSub), chan_sel_methods{iSel})
    end
end

for iSel=1:length(chan_sel_methods)
    plotSummaryGrand(chan_sel_methods{iSel})
end