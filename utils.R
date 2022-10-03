get_frex_report <- function(frex, frex_meterRel){
    isMeterRel = frex %in% frex_meterRel
    df2print <- data.frame(frex, isMeterRel, !isMeterRel)
    names(df2print) <- c('frequency (Hz)', 'is_meter_related', 'is_meter_unrelated')
    df2print$is_meter_related <- ifelse(df2print$is_meter_related==1, '✓', '')
    df2print$is_meter_unrelated <- ifelse(df2print$is_meter_unrelated==1, '✓', '')
    return(df2print)
}


prepare_frex <- function(df, frex, frex_meterRel, 
                         feat_to_zscore=c('amp_eeg', 'db_eeg', 'amp_coch', 'pow_coch', 'db_coch')){
    # keep only frequencies of interest 
    df <- df[df$freq %in% frex, ]
    # update which are meter-related frequencies 
    df$isMeterRel <- ifelse(df$freq %in% frex_meterRel, TRUE, FALSE)
    # get zscores 
    z_col_names <- sapply(feat_to_zscore, function(x) paste('z', x, sep='_'))
    df <- ddply(df, .(subject,rhythm,tone), function(d){
        for (i in c(1:length(feat_to_zscore))){
            col_name <- feat_to_zscore[i]
            new_col_name = z_col_names[i]
            d[, new_col_name] <- (d[, col_name] - mean(d[, col_name])) / sd(d[, col_name])
        }
        # d$z_eeg <- (d$amp_eeg-mean(d$amp_eeg))/sd(d$amp_eeg)
        # d$z_eeg_db <- (d$db_eeg-mean(d$db_eeg))/sd(d$db_eeg)
        # d$z_coch <- (d$amp_coch-mean(d$amp_coch))/sd(d$amp_coch)
        # d$z_coch_pow <- (d$pow_coch-mean(d$pow_coch))/sd(d$pow_coch)
        # d$z_coch_db <- (d$db_coch-mean(d$db_coch))/sd(d$db_coch)
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
    feat_str <- paste(z_col_names, collapse=' + ')
    f <- formula(paste(feat_str, '+ amp_eeg + amp_eeg_norm',
                       ' ~ ', 'subject + rhythm + tone + isMeterRel'))
    df_meter <- summaryBy(f, FUN=mean, data=df, keep.names=T)
    return(list(mean_amp=df_mean_amp, 
                meter=df_meter))
}


format_ttest_row <- function(res, d, ci){
    data.frame(t = res$statistic, 
               df = res$parameter, 
               mu = res$estimate,
               CI95_low = ci[1],
               CI95_high = ci[2],
               cohens_d = d$Cohens_d, 
               # d_CI_low = d$CI_low, 
               # d_CI_high = d$CI_high, 
               p = res$p.value, 
               signif = stars.pval(res$p.value))
}


boostrap_ci <- function(data, n=1000){
    b <- boot(data, 
              function(x, idx) mean(x[idx]), 
              n)
    ci <- boot::boot.ci(b, type='perc')
    return(ci$percent[4:5])
}


tVsCoch <- function(df, feature_eeg='z_amp_eeg', feature_coch='z_amp_coch'){
    res <- t.test(df[, feature_eeg], mu=df[1, feature_coch], alternative='greater')
    d <- cohens_d(df[, feature_eeg], mu=df[1, feature_coch], alternative='greater')
    ci <- boostrap_ci(df[, feature_eeg])
    format_ttest_row(res, d, ci)
}


tVs0 <- function(df, feature_eeg='z_amp_eeg'){
    res <- t.test(df[, feature_eeg], mu=0, alternative='greater')
    d <- cohens_d(df[, feature_eeg], mu=0, alternative='greater')
    ci <- boostrap_ci(df[, feature_eeg])
    format_ttest_row(res, d, ci)
    
}


tMeterRelvsUnrel <- function(df, feature='amp_eeg'){
    df_meter_rel <- filter(df, isMeterRel==TRUE) %>% arrange(subject)
    df_meter_unrel <- filter(df, isMeterRel==FALSE) %>% arrange(subject)
    stopifnot(length(unique(df_meter_rel$subject)) == nrow(df_meter_rel))
    stopifnot(length(unique(df_meter_unrel$subject)) == nrow(df_meter_unrel))
    val_meter_rel <- df_meter_rel[, feature]
    val_meter_unrel <- df_meter_unrel[, feature]
    res <- t.test(val_meter_rel, val_meter_unrel, paired=TRUE, alternative='greater')
    d <- cohens_d(val_meter_rel, val_meter_unrel, paired=TRUE, alternative='greater')
    ci <- boostrap_ci(val_meter_rel - val_meter_unrel)
    format_ttest_row(res, d, ci)
}


tLowVsHigh <- function(df, feature='amp_eeg'){
    df_low <- filter(df, tone=='low') %>% arrange(subject)
    df_high <- filter(df, tone=='high') %>% arrange(subject)
    stopifnot(length(unique(df_low$subject)) == nrow(df_low))
    stopifnot(length(unique(df_high$subject)) == nrow(df_high))
    val_low <- df_low[, feature]
    val_high <- df_high[, feature]
    res <- t.test(val_low, val_high, paired=TRUE, alternative='greater')
    d <- cohens_d(val_low, val_high, paired=TRUE, alternative='greater')
    ci <- boostrap_ci(val_low - val_high)
    format_ttest_row(res, d, ci)
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


get_meter_diff <- function(df, col_name_eeg){
    # calculate the difference meterRel - meterUnrel and check when it is positive vs. negative 
    # (we can use this to color the lines in a plot)
    get_diff <- function(df, feature){
        val_meter_rel <- filter(df, isMeterRel==TRUE) %>% select(.data[[feature]]) %>% unlist()
        val_meter_unrel <- filter(df, isMeterRel==FALSE) %>% select(.data[[feature]]) %>% unlist()
        val_diff <- val_meter_rel - val_meter_unrel
        data.frame(diff = val_diff,
                   diff_is_positive = ifelse(val_diff > 0, TRUE, FALSE))
    }
    df_diff <- ddply(df, 
                     .(subject, rhythm, tone), 
                     function(df) get_diff(df, feature=col_name_eeg)
    )
    df_diff$diff_is_positive <- factor(df_diff$diff_is_positive, 
                                       levels=c(TRUE, FALSE), 
                                       labels=c('rel > unrel', 'rel < unrel'))
    return(df_diff)
}


my_theme <- function(){
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
}


plot_z_meter <- function(df, col_name_eeg='z_amp_eeg', col_name_coch='z_amp_coch'){
    
    df$x <- as.numeric(df$tone)
    df$x_coch_start <- df$x-0.3
    df$x_coch_end <- df$x+0.3
    
    df_summary <- summarySEwithin(df, withinvars=c('rhythm','tone'), idvar='subject', measurevar=col_name_eeg)
    df_summary$x <- as.numeric(df_summary$tone)
    
    p <- ggplot(df, aes(x, .data[[col_name_eeg]], color=tone)) + 
        geom_segment(data=df[df$subject==df$subject[1],],
                     inherit.aes=FALSE,
                     aes(x=x_coch_start, xend=x_coch_end, 
                         y=.data[[col_name_coch]], yend=.data[[col_name_coch]],
                         color=tone),
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
        geom_errorbar(data=df_summary, 
                      aes(color=tone, ymin=.data[[col_name_eeg]]-ci, ymax=.data[[col_name_eeg]]+ci), 
                      size=1, width=0.2) + 
        scale_color_manual(name='tone', values=cols) + 
        facet_wrap(~rhythm) + 
        theme_cowplot() + 
        my_theme()
    
    return(p)
}


plot_amp_mean <- function(df, col_name_eeg='amp_eeg'){
    
    df$x <- as.numeric(df$tone)
    
    df_summary <- summarySEwithin(df, 
                                  withinvars=c('rhythm','tone'),
                                  idvar='subject', 
                                  measurevar=col_name_eeg)
    df_summary$x <- as.numeric(df_summary$tone)
    
    p <- ggplot(df, aes(x, .data[[col_name_eeg]], color=tone)) + 
        geom_line(col='grey80', aes(group=paste(subject,rhythm))) +
        geom_point(aes(group=paste(subject,rhythm)), size=2) +
        scale_color_manual(name='tone', values=cols_individual) + 
        scale_x_continuous(breaks=c(1:length(levels(df$tone))), labels=levels(df$tone), limits=c(1-0.3, 2+0.3)) +  
        new_scale_color() + new_scale_fill() + 
        geom_point(data=df_summary, aes(color=tone), size=3) + 
        geom_errorbar(data=df_summary,
                      aes(color=tone, 
                          ymin=.data[[col_name_eeg]]-ci, 
                          ymax=.data[[col_name_eeg]]+ci),
                      size=1, 
                      width=0.2) + 
        scale_color_manual(name='tone', values=cols) + 
        facet_wrap(~rhythm) + 
        theme_cowplot() + 
        my_theme()
    return(p) 
}


plot_amp_meter <- function(df, col_name_eeg='amp_eeg'){
    
    df$meter <- factor(df$isMeterRel, levels=c(TRUE, FALSE), labels=c('meter\nrel', 'meter\nunrel'))
    
    df_diff <- get_meter_diff(df, col_name_eeg)
    df <- merge(df, df_diff, by=c('subject', 'rhythm', 'tone'))
    
    p <- ggplot(df, aes(meter, .data[[col_name_eeg]])) + 
        geom_line(aes(color=diff_is_positive, group=paste(subject, rhythm, tone))) +
        scale_color_manual(name=NULL, values=c('#bdbbbb', '#666565')) + 
        new_scale_color() + new_scale_fill() + 
        geom_point(aes(color=tone, group=paste(subject, rhythm, tone)), size=2) +
        scale_color_manual(name='tone', values=cols_individual) + 
        facet_grid(rhythm ~ tone) + 
        theme_cowplot() + 
        my_theme() + 
        theme(
            axis.text.x = element_text(size=fontsize)
        )
    return(p) 
}


plot_amp_meter_diff <- function(df, col_name_eeg='amp_eeg'){
    df_diff <- get_meter_diff(df, col_name_eeg)
    # p <- ggplot(data=df_diff) + 
    #     geom_hline(yintercept=0, color='grey60', size=2) + 
    #     geom_half_violin(aes(x=rhythm, y=diff, split=tone, fill=tone), 
    #                      position='identity',
    #                      trim=FALSE) + 
    #     scale_fill_manual(name='tone', values=cols_individual) + 
    #     scale_x_discrete(name='rhythm') + 
    #     theme_cowplot() + 
    #     my_theme() + 
    #     theme(
    #         axis.text.x = element_text(size=fontsize),
    #         axis.title.x = element_text(size=fontsize)
    #     )
    p <- plot_amp_mean(df_diff, col_name_eeg='diff')
    return(p)
}

 
