# THIS IS SCRIPT FOR PLOTS
# the main script needs to be run before this script - as several dfs from the main script are needed

# LOG2 TRANSFORMATION PLOTS ----
# first raw plot of original values
raw_plot <- ggplot(df_long[which(!is.na(df_long$intensity)),], aes(x=intensity/10^6)) +
  geom_density(fill = "red", alpha = 0.65) +
  facet_wrap(~ sample) +
  scale_x_continuous(name = "quantity (*10^6)")+
  scale_fill_lancet()+
  theme(axis.text.x= element_text(angle = 90, vjust = 1, hjust=1 ,size =
                                    10)) +
  ggtitle("Density plot raw values")+
  gt_plot_theme

# then logged plot with log2 transformed values
log_plot <- ggplot(df_long[which(!is.na(df_long$intensity)),], aes(x=log2_intensity)) +
  geom_density(fill = "red", alpha = 0.65) +
  facet_wrap(~ sample) +
  scale_x_continuous(name = "quantity (log)")+
  theme(axis.text.x= element_text(angle = 90, vjust = 1, hjust=1 ,size =
                                    10)) +
  ggtitle("Density plot log2 values")+
  gt_plot_theme

# FILTERING PLOTS ----
# plot of NAs before filtering 
plot_before_filt <- ggplot(data = (df_long %>% group_by(protein) %>% summarise(n = sum(is.na(intensity)))),
                           aes(x = n))+
  geom_histogram(bins = 24, position = position_dodge(width = 0.5),
                  col = "dark blue", fill = "light blue", linewidth = 0.4)+
  ggtitle("Missing structure before filtering")+
  labs(x = "missing proteins",
       y = "samples (count)")+
  scale_x_continuous(breaks = c(0:23))+
  gt_plot_theme

# plot of NAs after filtering
plot_after_filt <- ggplot(data = (df_long_filt %>% group_by(protein) %>% summarise(n = sum(is.na(intensity)))),
                          aes(x = n))+
  geom_histogram(bins = 7, col = "dark blue", fill = "light blue")+
  ggtitle("Missing structure after filtering")+
  labs(x = "missing proteins",
       y = "samples (count)")+
  scale_x_continuous(breaks = c(0:6))+
  gt_plot_theme

# NORMALISATION PLOTS ----
# plot before normalisation
before_norm_plot <- ggplot(df_long_norm[which(!is.na(df_long_norm$intensity)),], aes(x = sample, y = log2_intensity))+
  geom_boxplot(aes(fill = group), alpha = 0.7)+
  ggtitle("Before normalisation")+
  scale_fill_jama()+
  labs(y = "Log2 intensity")+
  gt_plot_theme+
  theme(
    axis.text.x = element_text(angle = 270, hjust = 0))

# plot after normalisation
after_norm_plot <- ggplot(df_long_norm[which(!is.na(df_long_norm$intensity)),], aes(x = sample, y = norm_log2_intensity))+
  geom_boxplot(aes(fill = group), alpha = 0.7)+
  ggtitle("After normalisation")+
  scale_fill_jama()+
  labs(y = "Log2 intensity")+
  gt_plot_theme+
  theme(
    axis.text.x = element_text(angle = 270, hjust = 0))


# IMPUTATION plots----
# plot of valid results by sample and group after filtering but before imputation
valid_plot <- ggplot(df_long_filt %>% group_by(sample, group) %>% summarise(n = sum(!is.na(intensity))),
                     aes(x = sample, y = n, fill = group))+
  geom_col()+
  scale_fill_jama(alpha = 0.9)+
  ggtitle("Valid measurements after filtering per sample and group")+
  labs(y = "number of proteins")+
  gt_plot_theme+
  theme(axis.text.x = element_text(angle = 270, hjust = 0))


# plot displaying where imputation has been performed
impute_plot <- ggplot(df_long_imp, 
                      aes(x = norm_log2_intensity_imputed, fill = imputed, col = imputed))+
  geom_histogram(alpha = 0.85, bins = 60)+
  scale_color_manual(values = c("grey", "dark red"))+
  scale_fill_manual(values = c("grey", "dark red"))+
  facet_wrap(~sample)+
  ggtitle("Distribution of imputed values, by sample")+
  gt_plot_theme


# VOLCANO PLOT 
## volcano plot for LS vs Sepsis
# first determine diff expression for all proteins
preplot <- diff_expression(df_long_imp, protein_names = "protein", test_value = "norm_log2_intensity_imputed", 
                           group_var = "group", grp1 = "LS", grp2 = "Sepsis", sign_only = FALSE)

# then volcano plot
volcano_plot <- ggplot(preplot, aes(x = log2FC,y = -log10(qval), color = sign)) +
  geom_point(alpha = 0.9, size = 1) +
  geom_hline(yintercept = 1.3, linetype = 3, alpha = 0.5) +
  geom_vline(xintercept = 1.0, linetype = 3, alpha = 0.5) +
  geom_vline(xintercept = -1.0, linetype = 3, alpha = 0.5) +
  scale_x_continuous(limits = c(-4,4))+
  scale_y_continuous(limits = c(0,4))+
  scale_colour_manual(values = c("dark blue", "red")) +
  xlab("log2FC Other severe infections vs Lemierre´s syndrome") + ylab("-log10 q-value") +
  gt_plot_theme+
  theme(legend.position="none")


# HEATMAP ---
# heatmap preprocessing, obtaining protein names for significant
sign_protein_names <- data.frame(protein_id = table1[,"protein"],
                                 protein_name = sub("\\(.*", "", (GetNamesTaxa(table1[,"protein"]))$Protein.names)) %>%
  arrange(desc(protein_id))



# creating a preplot df
preplot <- df_long_imp %>% 
  filter(protein %in% sign_protein_names$protein_id,
         group %in% c("LS", "Sepsis")) %>% 
  group_by(protein, group) %>% 
  summarise(mean_intensity = mean(norm_log2_intensity_imputed)) 

# need to go via a wide dataframe to order proteins correctly
heat_plot_table <- pivot_wider(preplot, id_cols = protein, names_from=group, values_from=mean_intensity) %>%
  arrange(desc(protein))
heat_plot_table$name <- sign_protein_names$protein_name
heat_plot_table <- heat_plot_table %>% arrange(desc(LS))

heat_plot_table$name <- factor(heat_plot_table$name, levels = heat_plot_table$name)
# then to long df again
preplot <- pivot_longer(heat_plot_table, cols = 2:3, names_to = "group", values_to = "mean_intensity") %>% 
  mutate(group = ifelse(group == "Sepsis", "Other severe infections", "Lemierre´s syndrome"))

# this is the heatmap plot code
heat_plot <- ggplot(preplot, aes(x = group, y = fct_rev(name)))+
  geom_tile(aes(fill = mean_intensity))+
  scale_fill_viridis_c(option = "magma")+
  labs(y = NULL, x = NULL)+
  gt_plot_theme

heat_plot_table[,c(2,3)] <- round(heat_plot_table[,c(2,3)],2)
heat_plot_table <- heat_plot_table %>%
  rename(protein_id = protein, protein_name = name) %>%
  select(protein_id, protein_name, LS, Sepsis)


