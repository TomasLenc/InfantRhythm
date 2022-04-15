rm(list=ls())

reports_path <- '/datadisk/Dropbox/projects/XPInfant/derivatives/reports'

roi <- 'front'
rmarkdown::render('report_lmer.Rmd', output_file=file.path(reports_path, sprintf('roi-%s_reportLME.html',roi)))
rmarkdown::render('report_anova.Rmd', output_file=file.path(reports_path, sprintf('roi-%s_reportANOVA.html',roi)))

roi <- 'best'
rmarkdown::render('report_lmer.Rmd', output_file=file.path(reports_path, sprintf('roi-%s_reportLME.html',roi)))
rmarkdown::render('report_anova.Rmd', output_file=file.path(reports_path, sprintf('roi-%s_reportANOVA.html',roi)))
