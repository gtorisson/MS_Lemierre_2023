
# IMPUTATION
# this function imputes missing values as perseus default, with a random draw
# from a gaussian distribution
# input: df = a data frame in long format, sample = the variable (in quotes) designating unique samples,
# value = the value to be imputed, widht and downshift as Perseus defaults
impute_as_perseus <- function(df, sample, value, width = 0.3, downshift = 1.8){
  # need to create these before the loop
  df$imputed <- NA
  df$norm_log2_intensity_imputed <- NA
  # loop for each sample
  for(i in 1: length(unique(df[[sample]]))){
    # setting random seed for reproducability
    set.seed(124)
    # temp = temporary dataframe for each sample
    temp <- df %>% filter(as.numeric(sample)==i)
    temp_sd <- sd(temp[[value]], na.rm = T)
    temp_mean <- mean(temp[[value]], na.rm = T)
    # the sd and mean for imputation 
    small_sd <- width * temp_sd
    low_mean <- temp_mean - downshift * temp_sd
    # TRUE/FALSE to see where imputations are made
    temp$imputed <- is.na(temp[[value]])
    number_of_nas <- sum(temp$imputed)
    # first setting imputed to original values
    temp$norm_log2_intensity_imputed <- temp[[value]]
    # then imputing into where NA
    temp[which(temp$imputed),]$norm_log2_intensity_imputed <- 
      rnorm(number_of_nas, mean = low_mean, sd = small_sd)
    # incorporating into full long df for each sample
    df[which(as.numeric(df[[sample]])==i),]$imputed <- temp$imputed
    df[which(as.numeric(df[[sample]])==i),]$norm_log2_intensity_imputed <- temp$norm_log2_intensity_imputed
  }
  return(df)
}
 

# DIFF EXPRESSION
# this function creates a table of differentially expressed proteins
# df = a long data frame
# protein_names = the variables naming proteins
# test_value = the value to be tested
# group_var = the variable designating group allocation
# grp1 = the value of group_var designating the first group
# grp2 = the value of group_var designating the second group
# sign_only = if TRUE - only significantly differentially expressed proteins are returned

diff_expression <- function(df, protein_names, test_value, group_var, grp1, grp2,
                            sign_only = TRUE){
  #creating empty data frame with protein names
  results_df <- data.frame(protein = unique(df[[protein_names]]),
                           log2FC = NA, pval = NA)
  
  # this loop performs one-at-a-time students t-test for each protein between grp1 and grp2
  for(i in 1:nrow(results_df)){
    x = df[which(df[[group_var]] == grp1 & as.numeric(factor(df[[protein_names]])) == i),][[test_value]]
    y = df[which(df[[group_var]] == grp2 & as.numeric(factor(df[[protein_names]])) == i),][[test_value]]
    # returning p values and log2 fold change for each protein
   results_df$pval[i] <- (t.test(x, y, var.equal = T))$p.value
   results_df$log2FC[i] <- mean(x) - mean(y)
  }
  
  # further additions to results including q value and FC and significance
  results_df <- results_df %>% 
    mutate(qval = round(p.adjust(results_df$pval, method = "BH"),4),
           pval = as.numeric(format(round(pval,5), nsmall = 5)),
           FC = round(2^log2FC,2),
           sign = ifelse(((log2FC >= 1 |log2FC <=-1) & qval < 0.05), "+", "-"))
  
  if(sign_only){
    # only display significant results or not
    results_df <- results_df %>% filter(sign == "+")
  }
  # returning results arranged by decreasing log2FC
  return(results_df %>% arrange(desc(log2FC)))
}


# this is just a graphic theme for plots
gt_plot_theme <- theme_classic()+
  theme(axis.line = element_line(colour = "dark blue", linewidth = 0.25),
        axis.ticks = element_line(colour = "dark blue", linewidth = 0.25),
        axis.text = element_text(colour = "dark blue", size = 7),
        axis.title = element_text(colour = "dark blue", size = 9),
        legend.text = element_text(colour = "dark blue", size = 7),
        legend.title = element_text(colour = "dark blue", size = 8),
        title = element_text(colour = "dark blue", size = 11),
        strip.background = element_rect(fill = "light blue",
                                        linewidth = 0.75,color = "dark blue"),
        strip.text = element_text(size = 7, colour = "dark blue"))
