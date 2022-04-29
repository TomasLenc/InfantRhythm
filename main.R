rm(list=ls())

# ============== set this to your dataset root path ============== 
experiment_path <- '/datadisk/Dropbox/projects/XPInfant'
# ================================================================

reports_path <- file.path(experiment_path, 'derivatives/reports')
if (!dir.exists(reports_path)){
    dir.create(reports_path, recursive=TRUE)    
}

roi <- 'front'
rmarkdown::render('report_anova.Rmd', output_file=file.path(reports_path, sprintf('roi-%s_reportANOVA.html',roi)))

roi <- 'best'
rmarkdown::render('report_anova.Rmd', output_file=file.path(reports_path, sprintf('roi-%s_reportANOVA.html',roi)))
