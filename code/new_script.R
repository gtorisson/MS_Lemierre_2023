library(tidyverse) # data wrangling
library(readxl) # read in .xlsx files
library(stringi) # string manipulations
library(ggsci) # for plots in ggplot2
library(here) # cleaner file paths
library(UniprotR) # to retrieve protein names from Uniprot IDs

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

## setting the accepted fraction of missing per sample
missing_accepted = 0.3

# file import - this is relative to where the R project is
raw_data <- read_excel(here("data", "20220629_Karin_Holm_linkopingLib_Hong_Report_rawdata_nygruppering_avidentifierat.xlsx"), 2)
raw_data
names(raw_data)

# data is in wide format in the original data
df_wide <-  raw_data %>%
  # removing double names (separated by semicolons) in protein name
  mutate(PG.ProteinGroups = stri_extract_first_regex(PG.ProteinGroups, "^[^;]+")) %>%
  dplyr::rename(protein = PG.ProteinGroups) %>%
  # setting NA if filtered (low signal)
  mutate_all(~na_if(., "Filtered")) %>%
  # making readings from MS numeric (imported as character)
  mutate(across(2:30,~as.numeric(.))) %>%
  # setting NA if NaN 
  mutate_all(~ifelse(is.nan(.), NA, .))

# inspection - not necessary
df_wide

# making a long data frame 
df_long <- pivot_longer(df_wide, cols = c(2:30), names_to = c("group", "id"),
                        names_sep = -2, values_to = "intensity") %>% 
  # tidying the names from ... created by pivoting
  mutate_at(vars(2,3), ~str_remove_all(.,"[.]")) %>%
  # renaming for clarity
  mutate(group = ifelse(group == "NonLSthrombosis", "Other_thrombosis",group),
         # factorising
         group = factor(group),
         # naming individual samples
         sample = factor(str_c(group, id, sep = "_")),
         # creating a logged variable
         log2_intensity = log2(intensity)) %>% 
  # selecting only relevant variables
  select(protein, group,sample, intensity, log2_intensity)

# inspect
df_long

# plot of log2transformation 
# first raw plot
raw_plot <- ggplot(df_long[which(!is.na(df_long$intensity)),], aes(x=intensity/10^6)) +
  geom_density(fill = "red", alpha = 0.65) +
  facet_wrap(~ sample) +
  scale_x_continuous(name = "quantity (*10^6)")+
  scale_fill_lancet()+
  theme(axis.text.x= element_text(angle = 90, vjust = 1, hjust=1 ,size =
                                    10)) +
  ggtitle("Density plot raw values")+
  gt_plot_theme

raw_plot

# then logged plot
log_plot <- ggplot(df_long[which(!is.na(df_long$intensity)),], aes(x=log2_intensity)) +
  geom_density(fill = "red", alpha = 0.65) +
  facet_wrap(~ sample) +
  scale_x_continuous(name = "quantity (log)")+
  theme(axis.text.x= element_text(angle = 90, vjust = 1, hjust=1 ,size =
                                    10)) +
  ggtitle("Density plot log2 values")+
  gt_plot_theme

log_plot
# plot of NAs before filtering
plot_before_filt <- ggplot(data = (df_long %>% group_by(protein) %>% summarise(n = sum(is.na(intensity)))),
       aes(x = n))+
  geom_histogram( position = position_dodge(width = 0.5),
                  col = "dark blue", fill = "light blue", linewidth = 0.4)+
  ggtitle("Missing structure before filtering")+
  labs(x = "missing proteins",
       y = "samples (count)")+
  scale_x_continuous(breaks = c(0:29))+
  gt_plot_theme

plot_before_filt
## filtering away those with > 30% NA
df_long <- df_long %>% group_by(protein) %>% 
  filter(sum(is.na(intensity)) < missing_accepted*length(levels(df_long$sample)))

# plot of NAs after filtering
plot_after_filt <- ggplot(data = (df_long %>% group_by(protein) %>% summarise(n = sum(is.na(intensity)))),
       aes(x = n))+
  geom_histogram(bins = 9, col = "dark blue", fill = "light blue")+
  ggtitle("Missing structure after filtering")+
  labs(x = "missing proteins",
       y = "samples (count)")+
  scale_x_continuous(breaks = c(0:9))+
  gt_plot_theme

plot_after_filt
# normalisation using median subtraction per sample. Adding a normalised variable
df_long <- df_long %>% group_by (sample) %>% 
  mutate(norm_log2_intensity = log2_intensity - median(log2_intensity, na.rm = T))

df_long

# plot before normalisation
before_norm_plot <- ggplot(df_long[which(!is.na(df_long$intensity)),], aes(x = sample, y = log2_intensity))+
  geom_boxplot(aes(fill = group), alpha = 0.7)+
  ggtitle("Before normalisation")+
  scale_fill_jama()+
  labs(y = "Log2 intensity")+
  gt_plot_theme+
  theme(
    axis.text.x = element_text(angle = 270, hjust = 0))

before_norm_plot
# plot after normalisation - using the added variables
after_norm_plot <- ggplot(df_long[which(!is.na(df_long$intensity)),], aes(x = sample, y = norm_log2_intensity))+
  geom_boxplot(aes(fill = group), alpha = 0.7)+
  ggtitle("After normalisation")+
  scale_fill_jama()+
  labs(y = "Log2 intensity")+
  gt_plot_theme+
  theme(
    axis.text.x = element_text(angle = 270, hjust = 0))

after_norm_plot
# Imputation
# plot of valid results by sample and group
valid_plot <- ggplot(df_long %>% group_by(sample, group) %>% summarise(n = sum(!is.na(intensity))),
                     aes(x = sample, y = n, fill = group))+
  geom_col()+
  scale_fill_jama(alpha = 0.9)+
  ggtitle("Valid measurements after filtering per sample and group")+
  labs(y = "number of proteins")+
  gt_plot_theme+
  theme(axis.text.x = element_text(angle = 270, hjust = 0))

valid_plot
# start imputation by creating new variables and setting imputation parameters
df_long$imputed <- NA
df_long$norm_log2_intensity_imputed <- NA
width = 0.3
downshift = 1.8

# loop that for each sample (29) estimates the mean and sd of valid results and number of NAs, 
# then imputation with downshift and width and returning this to original df_long
for(i in 1 :length(unique(df_long$sample))){
  set.seed(124)
  temp <- df_long %>% filter(as.numeric(sample)==i)
  temp_sd <- sd(temp$norm_log2_intensity, na.rm = T)
  temp_mean <- mean(temp$norm_log2_intensity, na.rm = T)
  small_sd <- width * temp_sd
  low_mean <- temp_mean - downshift * temp_sd
  temp$imputed <- is.na(temp$norm_log2_intensity)
  number_of_nas <- sum(temp$imputed)
  temp$norm_log2_intensity_imputed <- temp$norm_log2_intensity
  temp[which(temp$imputed),]$norm_log2_intensity_imputed <- 
    rnorm(number_of_nas, mean = low_mean, sd = small_sd)
  df_long[which(as.numeric(df_long$sample)==i),]$imputed <- temp$imputed
  df_long[which(as.numeric(df_long$sample)==i),]$norm_log2_intensity_imputed <- temp$norm_log2_intensity_imputed
}

df_long
# plot displaying where imputation has been performed
impute_plot <- ggplot(df_long, 
       aes(x = norm_log2_intensity_imputed, fill = imputed, col = imputed))+
  geom_histogram(alpha = 0.85, bins = 60)+
  scale_color_manual(values = c("grey", "dark red"))+
  scale_fill_manual(values = c("grey", "dark red"))+
  facet_wrap(~sample)+
  ggtitle("Distribution of imputed values, by sample")+
  gt_plot_theme

impute_plot

## differential expression LS vs Sepsis
# creating an empty df
results_df <- data.frame(protein = unique(df_long$protein),
                         log2FC = NA, pval = NA, qval = NA)

# this loop performs a students t-test for each protein between LS and Sepsis
for(i in 1:nrow(results_df)){
  x = (df_long %>% filter(group == "LS" & as.numeric(factor(protein)) ==i))$norm_log2_intensity_imputed
  y = (df_long %>% filter(group == "Sepsis" & as.numeric(factor(protein)) ==i))$norm_log2_intensity_imputed
  results_df$pval[i] <- (t.test(x, y, var.equal = T))$p.value
  results_df$log2FC[i] <- mean(x) - mean(y)
}

# further additions to results including q value and FC and number of NAs
results_df <- results_df %>% 
  mutate(pval = as.numeric(format(round(pval,5), nsmall = 5)),
         qval = round(p.adjust(results_df$pval, method = "BH"),4),
         FC = round(2^log2FC,2),
         sign = ifelse(((log2FC >= 1 |log2FC <=-1) & qval < 0.05), "+", "-"),
         nas = (df_long %>% filter(group %in% c("LS", "Sepsis")) %>% 
                                     group_by(protein) %>% summarise(nas = sum(imputed)))$nas)

results_df

# results table, only significant are displayed (those with log2FC ≥±1 and q value < 0.05)
table1 <- results_df %>% filter(sign == "+") %>% 
  arrange(desc(log2FC))

# determining significant proteins and retrieving names
sign_proteins <- table1[,"protein"]
sign_proteins <- data.frame(protein = sign_proteins,
                       protein_name = sub("\\(.*", "", (GetNamesTaxa(sign_proteins))$Protein.names))


# results LS vs VTE - only estimated for proteins significant LS vs Sepsis
results_df2 <- data.frame(protein = unique(df_long$protein),
                         log2FC = NA, pval = NA, qval = NA)

# now this loop performs t-tests of LS vs DVT
for(i in 1:nrow(results_df2)){
  x = (df_long %>% filter(group == "LS" & as.numeric(factor(protein)) ==i))$norm_log2_intensity_imputed
  y = (df_long %>% filter(group == "DVT" & as.numeric(factor(protein)) ==i))$norm_log2_intensity_imputed
  results_df2$pval[i] <- (t.test(x, y, var.equal = T))$p.value
  results_df2$log2FC[i] <- mean(x) - mean(y)
}

# same here for LS vs DVT
results_df2 <- results_df2 %>% 
  mutate(pval = as.numeric(format(round(pval,5), nsmall = 5)),
         qval = round(p.adjust(results_df2$pval, method = "BH"),4),
         FC = round(2^log2FC,2),
         sign = ifelse(((log2FC >= 1 |log2FC <=-1) & qval < 0.05), "+", "-"),
         nas = (df_long %>% filter(group %in% c("LS", "DVT")) %>% 
                  group_by(protein) %>% summarise(nas = sum(imputed)))$nas)

table2 <- results_df2 %>% filter(protein %in% sign_proteins$protein)

table2
# combinina into a new table with both comparisons
table3 <- left_join((table1 %>% select(protein, FC, qval, sign) %>% rename(FC1 = FC, q1 = qval, sign1 = sign)),
          (table2 %>% select(protein, FC, qval, sign) %>% rename(FC2 = FC, q2 = qval, sign2 = sign)))


table3
## volcano plot for LS vs Sepsis
volcano_plot <- ggplot(results_df, aes(x = log2FC,y = -log10(qval), color = sign)) +
  geom_point(alpha = 0.9, size = 1) +
  geom_hline(yintercept = 1.3, linetype = 3, alpha = 0.5) +
  geom_vline(xintercept = 1.0, linetype = 3, alpha = 0.5) +
  geom_vline(xintercept = -1.0, linetype = 3, alpha = 0.5) +
  scale_x_continuous(limits = c(-4,4))+
  scale_y_continuous(limits = c(0,4))+
  scale_colour_manual(values = c("dark blue", "red")) +
  xlab("log2FC Other severe infections vs Lemierre") + ylab("-log10 q-value") +
  gt_plot_theme+
  theme(legend.position="none")

volcano_plot
## heatmap preprocessing
preplot <- df_long %>% 
  filter(protein %in% sign_proteins$protein,
         group %in% c("LS", "Sepsis")) %>% 
           group_by(protein, group) %>% 
  summarise(mean_intensity = mean(norm_log2_intensity_imputed))

# need to go via a wide dataframe to order proteins correctly
temp <- pivot_wider(preplot, id_cols = protein, names_from=group, values_from=mean_intensity)
temp <- temp %>% arrange(desc(LS))
temp$protein <- factor(temp$protein, levels = temp$protein)
temp$name <- sign_proteins$protein_name
temp$name <- factor(temp$name, levels = temp$name)
# then long file again
preplot <- pivot_longer(temp, cols = 2:3, names_to = "group", values_to = "mean_intensity")


# this is the heatmap code
heat_plot <- ggplot(preplot, aes(x = group, y = fct_rev(name)))+
  geom_tile(aes(fill = mean_intensity))+
  scale_fill_viridis_c(option = "magma")+
  labs(y = NULL)+
  gt_plot_theme


heat_plot
