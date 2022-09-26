get_frex_report <- function(frex, frex_meterRel){
    isMeterRel = frex %in% frex_meterRel
    df2print <- data.frame(frex, isMeterRel, !isMeterRel)
    names(df2print) <- c('frequency (Hz)', 'is_meter_related', 'is_meter_unrelated')
    df2print$is_meter_related <- ifelse(df2print$is_meter_related==1, '✓', '')
    df2print$is_meter_unrelated <- ifelse(df2print$is_meter_unrelated==1, '✓', '')
    return(df2print)
}

prepare_frex <- function(df, frex, frex_meterRel){
    # prepare mean meter-related frequencies 
    which_meterRel <- which(frex %in% frex_meterRel)
    # keep only frequencies of interest 
    df <- df[df$freq %in% frex, ]
    # update which are meter-related frequencies 
    df$isMeterRel <- ifelse(df$freq %in% frex_meterRel, TRUE, FALSE)
    # get zscores 
    df <- ddply(df, .(subject,rhythm,tone), function(d){
        d$z_eeg <- (d$amp_eeg-mean(d$amp_eeg))/sd(d$amp_eeg)
        d$z_eeg_db <- (d$db_eeg-mean(d$db_eeg))/sd(d$db_eeg)
        d$z_coch <- (d$amp_coch-mean(d$amp_coch))/sd(d$amp_coch)
        d$z_coch_pow <- (d$pow_coch-mean(d$pow_coch))/sd(d$pow_coch)
        d$z_coch_db <- (d$db_coch-mean(d$db_coch))/sd(d$db_coch)
        return(d)
    })
    # get normalized amplitudes (from 0 to 1)
    df <- ddply(df, .(subject,rhythm,tone), function(d){
        d$amp_eeg_norm <- (d$amp_eeg-min(d$amp_eeg)) / max((d$amp_eeg-min(d$amp_eeg)))
        return(d)
    })    
    stopifnot(min(df$amp_eeg_norm)==0)
    stopifnot(max(df$amp_eeg_norm)==1)
    # take mean amplitude across all frequencies of interest 
    df_mean_amp <- summaryBy(amp_eeg ~ subject+rhythm+tone, FUN=mean, data=df, keep.names=T)
    # take mean zscore across meter-related frequencies 
    df_meter <- summaryBy(amp_eeg + amp_eeg_norm + z_eeg + z_eeg_db + z_coch + z_coch_pow + z_coch_db ~
                              subject + rhythm + tone + isMeterRel, 
                          FUN=mean, data=df, keep.names=T)
    return(list(mean_amp=df_mean_amp, 
                meter=df_meter))
}

tVsCoch <- function(df, feature_eeg='z_eeg', feature_coch='z_coch'){
    res <- t.test(df[, feature_eeg], mu=df[1, feature_coch], alternative='greater')
    d <- cohens_d(df[, feature_eeg], mu=df[1, feature_coch], alternative='greater')
    data.frame(t = res$statistic, 
               df = res$parameter, 
               cohens_d = d$Cohens_d, 
               d_CI_low = d$CI_low, 
               d_CI_high = d$CI_high, 
               p = res$p.value, 
               signif = stars.pval(res$p.value))
}

tVs0 <- function(df){
    res <- t.test(df$z_eeg, mu=0, alternative='greater')
    d <- cohens_d(df$z_eeg, mu=0, alternative='greater')
    data.frame(t = res$statistic, 
               df = res$parameter, 
               cohens_d = d$Cohens_d, 
               d_CI_low = d$CI_low, 
               d_CI_high = d$CI_high, 
               p = res$p.value, 
               signif = stars.pval(res$p.value))
}

tMeterRelvsUnrel <- function(df, feature='amp_eeg'){
    val_meterRel <- df[df$isMeterRel==TRUE,feature]
    val_meterUnrel <- df[df$isMeterRel==FALSE,feature]
    res <- t.test(val_meterRel, val_meterUnrel, paired=TRUE, alternative='greater')
    d <- cohens_d(val_meterRel, val_meterUnrel, paired=TRUE, alternative='greater')
    data.frame(t = res$statistic, 
               df = res$parameter, 
               d = d$Cohens_d, 
               d_CI_low = d$CI_low, 
               d_CI_high = d$CI_high, 
               p = res$p.value, 
               signif = stars.pval(res$p.value))
}

format_anova_table <- function(model){
    res <- model$ANOVA
    res$F <- round(res$F,2)
    res$p <- round(res$p,5)
    res$ges <- NULL
    
    esq <- eta_squared(model$aov, partial=TRUE)
    res$`η2(partial)` <- round(esq$Eta2_partial, 2)
    res$`η2_CI_low` <- round(esq$CI_low, 2)
    res$`η2_CI_high` <- round(esq$CI_high, 2)
    
    res <- res %>% 
        relocate(c('p', 'p<.05'), .after=last_col())
    
    return(res)
}

get_pairwise_contrasts <- function(model){
    # calcualte contrast
    emm <- emmeans(model, 'tone')
    c <- contrast(emm, 'pairwise', adjust='none')
    # get pvalues
    c_pval <- as.data.frame(c)
    c_pval$p.value <- p.adjust(c_pval$p.value, adjust_method)
    c_pval$significance <- stars.pval(c_pval$p.value)
    options(scipen=100)
    c_pval$p.value <- format.pval(c_pval$p.value,digits=1,eps=0.0001)
    # get CIs
    c_ci <- as.data.frame(confint(c, adjust=adjust_method))
    # merge it together
    c_merged <- merge(c_pval, c_ci %>% dplyr::select(contrast,lower.CL,upper.CL), by=c('contrast'))
    c_merged <- c_merged %>% 
        dplyr::select(contrast,estimate,df,t.ratio,lower.CL,upper.CL,p.value,significance)
    c_merged <- c_merged %>% 
        mutate_at(c('estimate','t.ratio','lower.CL','upper.CL'), round, digits=3)
    return(c_merged)
}

plot_z_meter <- function(df){
    
    df$x <- as.numeric(df$tone)
    df$x_coch_start <- df$x-0.3
    df$x_coch_end <- df$x+0.3
    
    df_summary <- summarySEwithin(df, withinvars=c('rhythm','tone'), idvar='subject', measurevar='z_eeg')
    df_summary$x <- as.numeric(df_summary$tone)
    
    p <- ggplot(df, aes(x, z_eeg, color=tone)) + 
        geom_segment(data=df[df$subject==df$subject[1],],
                     inherit.aes=FALSE,
                     aes(x=x_coch_start, xend=x_coch_end, y=z_coch, yend=z_coch, color=tone),
                     size=3, alpha=0.5) +
        scale_color_manual(name='tone', values=cols) + 
        new_scale_color() + new_scale_fill() + 
        geom_line(col='grey80', aes(group=paste(subject,rhythm))) +
        geom_point(aes(color=tone), size=2) +
        scale_color_manual(name='tone', values=cols_individual) + 
        scale_x_continuous(breaks=c(1:length(levels(df$tone))), labels=levels(df$tone)) +  
        geom_hline(yintercept=0) +
        new_scale_color() + new_scale_fill() + 
        geom_point(data=df_summary, aes(color=tone), size=3) + 
        geom_errorbar(data=df_summary, aes(color=tone, ymin=z_eeg-ci, ymax=z_eeg+ci), size=1, width=0.2) + 
        scale_color_manual(name='tone', values=cols) + 
        facet_wrap(~rhythm) + 
        theme_cowplot() + 
        theme(
            axis.line.x = element_blank(), 
            axis.text.x = element_blank(), 
            axis.text.y = element_text(size=fontsize), 
            axis.title = element_text(size=fontsize), 
            axis.ticks = element_blank(), 
            axis.title.x = element_blank(),
            legend.text = element_text(size=fontsize), 
            legend.title = element_text(size=fontsize), 
            strip.background = element_blank(), 
            strip.text = element_text(size=fontsize, face='bold')
        )
    return(p)
}

plot_amp_mean <- function(df){
    
    df$x <- as.numeric(df$tone)
    
    df_summary <- summarySEwithin(df, withinvars=c('rhythm','tone'), idvar='subject', measurevar='amp_eeg')
    df_summary$x <- as.numeric(df_summary$tone)
    
    p <- ggplot(df, aes(x, amp_eeg, color=tone)) + 
        geom_line(col='grey80', aes(group=paste(subject,rhythm))) +
        geom_point(aes(group=paste(subject,rhythm)), size=2) +
        scale_color_manual(name='tone', values=cols_individual) + 
        scale_x_continuous(breaks=c(1:length(levels(df$tone))), labels=levels(df$tone), limits=c(1-0.3, 2+0.3)) +  
        new_scale_color() + new_scale_fill() + 
        geom_point(data=df_summary, aes(color=tone), size=3) + 
        geom_errorbar(data=df_summary, aes(color=tone, ymin=amp_eeg-ci, ymax=amp_eeg+ci), size=1, width=0.2) + 
        scale_color_manual(name='tone', values=cols) + 
        facet_wrap(~rhythm) + 
        theme_cowplot() + 
        theme(
            axis.text = element_text(size=fontsize), 
            axis.title = element_text(size=fontsize), 
            axis.ticks = element_blank(), 
            axis.title.x = element_blank(),
            legend.text = element_text(size=fontsize), 
            legend.title = element_text(size=fontsize), 
            strip.background = element_blank(), 
            strip.text = element_text(size=fontsize, face='bold')
        )
    return(p) 
}












 
