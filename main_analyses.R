###############################
###############################
#####  SCALING - WISC FSIQ ####
###############################
###############################



##############
### SET UP ###
##############

library(tidyverse)
library(magrittr)
library(ComplexHeatmap)
library(reshape)
library(reshape2)
library(ggpmisc)
library(ggpubr)
library(visdat)
library(naniar)
library(broom)
library(tableone)
library(finalfit)
library(splines)
library(mgcv)
library(psych)
library(RColorBrewer)
# use theme_bw() for all plots
theme_set(theme_bw())
setwd("~/Desktop/HBN") # path to project folder



#############################
### READ IN CLEANED  DATA ###
#############################

subscales <- read.csv("./modified_data/ninth_release_iq_analyses2.csv")



##################################
### SPLITTING IQ INTO TERTILES ###
##################################

# visualize distribution
ggplot(subscales, aes(WISC_FSIQ)) + geom_histogram() + labs(x = "IQ Score", y = "Count")

# dfs of top and bottom 33% of IQ distribution
top33 <- subscales %>%
  filter(!is.na(WISC_FSIQ)) %>%    # remove obs with no WISC IQ data
  slice_max(., WISC_FSIQ, prop = 0.33)
bottom33 <- subscales %>%
  filter(!is.na(WISC_FSIQ)) %>%    # remove obs with no WISC IQ data
  slice_min(., WISC_FSIQ, prop = 0.33)
# collapse into a list
split <- list(high = top33, low = bottom33)

# get an idea of how many observations are NA for each subscale in each tertile
compare <- data.frame(high = apply(top33, 2, function(x) sum(!is.na(x))), low = apply(bottom33, 2, function(x) sum(!is.na(x))))
### all have at least 300 obs in each tertile



############################
### IQ X SUBSCALE PLOTS  ###  
############################

# pivot to long form for plotting
subscales.long <- subscales %>%
  # only take the subscales, not the extra demographic info
  select(CBCL_AD_T : WISC_FSIQ) %>%
  pivot_longer(., cols = !WISC_FSIQ, names_to = "subscale", values_to = "score") %>%
  filter(str_detect(subscale, "WISC") == FALSE)

# df with just the iq x subscale correlations, to use in plots
cor_with_iq <- data.frame(cor(subscales[17:92], use = "pairwise.complete.obs")) %>%
  select(WISC_FSIQ) %>% 
  dplyr::rename(corr = WISC_FSIQ) %>%
  rownames_to_column(var = "subscale")


# save this df for later comparison with other scaling vars
write.csv(cor_with_iq, file = "./modified_data/corr_with_WISC_FSIQ_33.csv")

# function for plotting
PLOT_FXN <- function(x) {
  # plot wisc against subscale score
  output <- ggplot(x, aes(WISC_FSIQ, score)) +
    geom_point(aes(col = subscale), show.legend = FALSE, alpha = 0.3) +
    # add regression line for each subscale
    geom_smooth(method = "lm", aes(color = subscale), show.legend = FALSE) +
    # separate facet for each subscale
    facet_wrap(~ subscale) +
    # make facet labels slightly larger
    theme(strip.text = element_text(size = 13)) +
    # add r correlation
    geom_text(data = cor_with_iq[cor_with_iq$subscale %in% x[["subscale"]], , ], aes(x = Inf, y = Inf, hjust = 1, vjust = 1.5,
                                                                                     label = paste0("r = ", round(corr, 4))), size = 4) +
    # add annotation printing the p-value for each relationship
    stat_fit_glance(method = "lm",
                    method.args = list(formula = y ~ x), 
                    label.y = "top", 
                    mapping = aes(label = sprintf('~italic(P)~"="~%.2g',
                                                  stat(p.value))), 
                    parse = TRUE)
  return(output)
}

# apply this function to groups of subscales to plot
( plts <- subscales.long %>%
    # split into separate df for each instrument (can be delineated based on beginning of colname), all contained in a list
    group_split(substr(subscales.long$subscale, 1, 4)) %>%
    # then map the PLOT_FXN onto each df in this list
    map(., PLOT_FXN) )

# p-values of correlations
corr_p_val <- apply(subscales[, 17:86], 2, FUN = function(c) cor.test(x = c, y = subscales$WISC_FSIQ)) %>%
  map_dfc(., function(d){
    d[["p.value"]]
  }) %>%
  t() %>%
  data.frame(p_val = .) %>%
  rownames_to_column(var = "subscale")

# put r and p values into same df
assertthat::assert_that(identical(corr_p_val$subscale, cor_with_iq$subscale[1:70]))
iq_cor_export <- cor_with_iq %>%
  magrittr::extract(1:70, ) %>%
  mutate(corr = round(corr, 3),
         p_val = round(corr_p_val$p_val, 3))
# save 
write.csv(iq_cor_export, file = "./modified_data/iq_corr_p_val_fortable_33.csv")

# sig at alpha = 0.05; 47
corr_p_val %>%
  filter(p_val <= 0.05)
# sig after Bonferroni; 32
corr_p_val %>%
  filter(p_val <= (0.05/70))


#####################################################
### MEAN SUBSCALE SCORES IN TOP VS BOTTOM TERTILE ###
#####################################################

# calculate mean score for each subscale in top33 and bottom33 
observed_means <- map_dfc(split, function(d){
  d %>%
    # select subscales
    select(CBCL_AD_T : ICU_SR_Tot) %>%
    # compute mean
    apply(., 2, mean, na.rm = TRUE) %>%
    data.frame(.) %>%
    magrittr::set_rownames(colnames(d[17:86]))
}) %>%
# save the means data frame in this order, so i can use it to compare against the null distribution later on
  set_colnames(c("high", "low")) %>%
  rownames_to_column("subscale") %>%
  # add a var for the difference between mean score in low vs high IQ groups
  mutate(difference = high - low) 



#########################################################
### COMPARING MEAN DIFFERENCES TO A NULL DISTRIBUTION ###
#########################################################

NULL_MEANS <- function(x){   # input x is a df in which rows are observations, columns are subscales
  # randomly generate vector of row numbers that is double the length of the tertile splits
  sample_size <- floor(0.66*nrow(x))
  selected <- sample(x = nrow(x), size = sample_size)
  # then take the top and bottom half of this random vector of indices
  inda <- head(selected, n = (0.5 * length(selected)))
  indb <- tail(selected, n = (0.5 * length(selected)))
  # split df based on these indices, so now we have two dfs of a random subset of 30% of the observations
  terta <- x[inda, ]
  tertb <- x[indb, ]
  # calculate the means for each subscale
  meansa <- apply(terta, 2, mean, na.rm = TRUE)
  meansb <- apply(tertb, 2, mean, na.rm = TRUE)
  # verify row and column names are in the same order
  assertthat::assert_that(identical(names(meansa), names(meansb)))
  # subtract one matrix from the other to get difference matrix
  diff <- meansa - meansb
  return(diff)
}
# iterate the null function 1000 times over the "subscales" dataframe, to get the null distribution of the difference in means for each subscale
n <- 10000
null_means <- map(seq_len(n), ~ NULL_MEANS(x = subscales[, 17:86])) 
# give each df a name other than a number, so that I can use the df names as column names when I collapse the list
names(null_means) <- paste("V", seq_len(n), sep = "")
# make df in which rows are subscales and columns are each iteration of the null means fxn
null_means.vec <- Reduce(cbind, null_means) %>%
  data.frame(.) %>%
  rownames_to_column(.data = ., var = "subscale")


## FIND THE PERCENTILE AT WHICH EACH OBSERVED VALUE FALLS ##
# check that subscales are in the same order in observed and null data sets
assertthat::assert_that(identical(observed_means$subscale, null_means.vec$subscale))

for (i in 1:nrow(observed_means)){
  observed_means$percentile[i] <- ecdf(null_means.vec[i, 2:10001])(observed_means$difference[i])
}
# NOTE: percentile values will be slightly different with each time the code is run, because the randomized row groupings
# change each time; however, the observed value will not change

# display table of significant subscales
observed_means %>%
  filter(percentile >= 0.975 | percentile <= 0.025) %>% # sig at alpha = 0.05
  arrange(desc(percentile))  # arrange in order of decreasing significance

# how many sig after Bonferroni corrected
observed_means %>% 
  filter(percentile <= 0.0003571428 | percentile >= 0.9996429)

# create var for p-value in addition to percentile
observed_means <- observed_means %>%
  mutate(p_value = ifelse(test = percentile > 0.5, yes = 1 - percentile, no = percentile), 
         bf_p_value = ifelse(test = (p_value * 70) <= 1, yes = p_value * 70, no = 1))

# save observed_means for comparison
write.csv(observed_means, file = "./modified_data/diff_in_means_WISC_FSIQ_33.csv")
# save the null distribution of differences in means to use for all the other WISC subscales
write.csv(null_means.vec, file = "./modified_data/null_means_distribution_33.csv")

### BAR PLOT SHOWING DIFFERENCE IN MEANS ###
observed_means %>%
  mutate(sig = ifelse(test = percentile >= 0.975 | percentile <= 0.025, yes = TRUE, no = FALSE)) %>%
  ggplot(., aes(fct_reorder(subscale, desc(difference)), difference)) +
  geom_col(aes(fill = sig)) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "Subscale", y = "Difference in mean score (low - high)", fill = "Significant")



#####################################
### RANKING COHEN'S D EFFECT SIZE ###
#####################################

# calculate standard deviation of each subscale in the high and low IQ groups
stdev <- map_dfc(split, function(d){
  d %>%
    # select subscales
    select(CBCL_AD_T : ICU_SR_Tot) %>%
    # compute standard deviation
    apply(., 2, sd, na.rm = TRUE) %>%
    # make df
    data.frame(.) %>%
    # keep subscales as row names
    magrittr::set_rownames(colnames(d[17:86]))
}) %>%
  # edits to colnames and save subscales as a column, not rownames
  set_colnames(c("high", "low")) %>%
  rownames_to_column("subscale") %>%
  # compute the pooled standard deviation, which will be used to caluclate cohen's d
  mutate(pooled_sd = sqrt(((high ^ 2) + (low ^ 2)) / 2) )

# divide the difference in means by the pooled standard deviation
# first check that rows are in the same order
assertthat::assert_that(identical(observed_means$subscale, stdev$subscale))
cohens_d <- data.frame(subscale = observed_means$subscale, d = observed_means$difference / stdev$pooled_sd) %>%
  mutate(sig = ifelse(observed_means$percentile > 0.975 | observed_means$percentile < 0.025, # if significant at alpha = 0.05
                      TRUE, FALSE))

# present in order of highest to lowest effect size
cohens_d %>%
  arrange(desc(d)) 

# save cohen's d so that I can compare to splits with other cogx scaling vars
write.csv(cohens_d, file = "./modified_data/cohens_d_WISC_FSIQ_33.csv")

which(observed_means$p_value <= (0.025/70)) %>% length  # 29 survive Bonferroni


### BAR PLOT SHOWING COHEN'S D ###
# more human-readable names for variables, for figures in paper
new_names <- c("Anxious Dep", "Withdrawn Dep", "Somatic", "Social", "Thought", 
               "Attention", "Rule Breaking", "Aggression", "Internalizing", "Externalizing", 
               "Total", "Aggression", "Anxious Dep", "Attention", "Withdrawn Dep", 
               "Rule Breaking", "Somatic", "Social", "Thought", "Externalizing", "Internalizing", 
               "Total", "Aggression", "Relationships", "Hyperactivity", "Inattentiveness", "Learning", 
               "Social Awareness", "Social Cognition", "Social Communication", "Repetitive Behaviors", 
               "Social Motivation", "Repetitive Behaviors", "Social Comm/Interaction", "Total Social",
               "Generalized Anx", "Panic", "Social Anx", "School Avoidance", "Separation Anx", "Total Anx", 
               "Generalized Anx", "Panic", "Social Anx", "School Avoidance", "Separation Anx", "Total Anx",
               "Stereotyped Behvaiors", "Self-injury", "Compulsive Behaviors", "Ritualized Behaviors", 
               "Restricted Interests", "Total Repetitive", "Conduct", "Total", "Emotional", "Externalizing",
               "Total", "Hyperactive", "Internalizing", "Peer Problems", "Prosocial (reverse)",
               "Callousness", "Uncaring", "Unemotional", "Total Psychopathy", "Callousness", 
               "Uncaring", "Unemotional", "Total Psychopathy")
current_names <- as.vector(observed_means$subscale)
# adding column for instrument labels
observed_means <- observed_means %>%
  mutate(instrument = c(rep("CBCL", 11), rep("YSR", 11), rep("CSR", 5), rep("SRS", 8), rep("SCARED_PR", 6), 
                        rep("SCARED_SR", 6), rep("RBS", 6), rep("SDQ", 9), rep("ICU_PR", 4), rep("ICU_SR",4)), 
         cohens_d = cohens_d$d)
# instrument color assignments
instrument_colors <- structure(brewer.pal(n = 10, name = "Set3"), 
                               names = c("CBCL", "YSR", "CSR", "SRS", "SCARED_PR", "SCARED_SR", "RBS", "SDQ", "ICU_PR", "ICU_SR")) 
inst_colors_df <- observed_means %>%
  select(subscale, instrument) %>%
  mutate(color = instrument_colors[instrument], 
         hr_name = new_names)
inst_colors_vec <- structure(inst_colors_df$color, names = inst_colors_df$subscale)
# save
write.csv(instrument_colors, file = "./modified_data/instrument_color_assignments10.csv")
write.csv(inst_colors_df, file = "./modified_data/instrument_color_assignments_df.csv")
write.csv(inst_colors_vec, file = "./modified_data/instrument_color_assignments70.csv")
# cohen's d plot
( cohens_d_plot <- observed_means %>%
  # change sig for Bonferroni corrected p values
  mutate(sig = ifelse(test = percentile >= 0.9996429 | percentile <= 0.0003571428, yes = TRUE, no = FALSE)) %>%
  ggplot(., aes(fct_reorder(subscale, cohens_d), cohens_d)) +
  geom_col(aes(fill = sig)) +
  theme(axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5), 
        text = element_text(size = 13.5)) +
  labs(x = "Subscale", y = "Cohen's *d*", fill = "Significant") +
  # change colors so they're not the same as tertile split colors
  scale_fill_manual(values=c("#999999", "#56B4E9")) +
  # put the subscale names on the top so it's a little easier to interpret
  scale_x_discrete(position = "top", 
                   # replace variable names with human-readable names
                   breaks = current_names, labels = new_names) +
  # need this command for making d italic in y axis
  theme(axis.title.y = ggtext::element_markdown(), 
        # put legend inside of plot area
        legend.position = c(0.8, 0.2)) )

# adding annotation bar for instrument 
legend_x <- ggplot(observed_means, aes(x = fct_reorder(subscale, cohens_d), y = 0.1)) + 
  geom_point(aes(color = instrument), shape = 15, size = 3, show.legend = FALSE) + 
  theme_classic() + 
  theme(axis.title = element_blank(), axis.line = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), 
        plot.margin = unit(c(0,0,0,0), "cm")) +
  labs(color = "Instrument Name") +
  scale_color_manual(values = instrument_colors)
cohens_d_plot + annotation_custom(ggplotGrob(legend_x),
                                  ymin = 0.1, ymax = 0.13, 
                                  xmin = 0, xmax = 70.7) 


######################################
### TERTILE-WISE DIFFERENCE MATRIX ###
######################################

# calculate correlation matrix for each tertile, then z-transform
split.corr <- map(split, function(d){
  d %>% 
    select(CBCL_AD_T : ICU_SR_Tot) %>%
    cor(., use = "pairwise.complete.obs") %>% 
    # z-transform
    fisherz()
})

# subtract LOW from HIGH to make difference matrix, in which positive difference means that pair of subscales are more highly
# correlated in HIGH IQ than LOW IQ individuals -- NOTE this is the reverse of the difference matrices I made previously, 
# but I'm making them this way now so the colors are consonant with the dot-product method matrix
assertthat::assert_that(identical(colnames(split.corr[["low"]]),  colnames(split.corr[["high"]]))) # check row and col names are identical
assertthat::assert_that(identical(rownames(split.corr[["low"]]),  rownames(split.corr[["high"]])))
# make diagonal values 0 instead of Inf
split.corr[[1]][which(is.infinite(split.corr[[1]]) == TRUE)]<- 0
split.corr[[2]][which(is.infinite(split.corr[[2]]) == TRUE)]<- 0
# subtract
lmh_z <- split.corr[["high"]] - split.corr[["low"]]  
# heatmap
fixedcolors <- circlize::colorRamp2(breaks = c(0.3, 0, -0.3), colors = c("#276419", "white", "#8E0152"))
col2 <- circlize::colorRamp2(breaks = c(1, 0, -1), colors = c("#276419", "white", "#8E0152"))
# assign rows and columns to a certain order (determined by hierarchical clustering) so that I can match it in the NS matrix
clusters <- hclust(d = dist(x = lmh_z, method = "euclidean"), method = "complete")
lmh_reordered <- data.frame(lmh_z) %>%
  select(clusters[["order"]]) %>%
  t() %>%
  data.frame() %>%
  select(clusters[["order"]]) %>%
  as.matrix()
# human-readable names for figure
new_names_heatmap <- structure(new_names, names = colnames(lmh_z))
# visualize
( hm_diff <- Heatmap(lmh_reordered, name = "Tertile difference matrix", 
        col = fixedcolors, 
        row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 8),
        # human-readable names
        row_labels = new_names_heatmap[colnames(lmh_reordered)], column_labels = new_names_heatmap[colnames(lmh_reordered)],
        # square
        height = unit(16.5, "cm"), width = unit(16.5, "cm"),
        cluster_rows = FALSE, cluster_columns = FALSE) )
# distribution of differences in pairwise correlation
hist(lmh_z)
# avg of all edges
mean(lmh_z[upper.tri(lmh_z, diag = FALSE)])

# save to compare
write.csv(lmh_z, file = "./modified_data/diff_mat_WISC_FSIQ_33_ztransform.csv")
write.csv(new_names_heatmap, file = "./modified_data/humanreadablenames.csv")

# now with instrument colors
inst_legend1 <- Heatmap(as.matrix(colnames(lmh_reordered)), name = "Instrument Name", 
                        col = inst_colors_vec, cluster_rows = FALSE, 
                        show_heatmap_legend = FALSE, height = unit(16.5, "cm"), width = unit(0.5, "cm"))
fig_diffmat <- hm_diff + inst_legend1 + 
  # add in names bc adding 2nd heatmap deletes them
  rowAnnotation(rn = anno_text(new_names_heatmap[colnames(lmh_reordered)]))
draw(fig_diffmat)


### CORRELATION MATRICES IN EACH TERTILE INDIVIDUALLY ###
# assign rows and columns to a certain order (determined by hierarchical clustering) so that I can add color annotation for instruments
clusters_high <- hclust(d = dist(x = split.corr[[1]], method = "euclidean"), method = "complete")
upper_tert_reordered <- data.frame(split.corr[[1]]) %>%
  select(clusters_high[["order"]]) %>%
  t() %>%
  data.frame() %>%
  select(clusters_high[["order"]]) %>%
  as.matrix()
( hm_upper <- Heatmap(upper_tert_reordered, name = "Upper tertile z matrix", 
        col = col2, show_row_dend = FALSE, show_column_dend = FALSE,
        cluster_rows = FALSE, cluster_columns = FALSE,
        # human-readable names
        row_labels = new_names_heatmap[colnames(upper_tert_reordered)], column_labels = new_names_heatmap[colnames(upper_tert_reordered)],
        # square
        height = unit(16.5, "cm"), width = unit(16.5, "cm"),
        row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 7)) )
# now with instrument colors
inst_legend_upper <- Heatmap(as.matrix(colnames(upper_tert_reordered)), name = "Instrument Name",
                        col = inst_colors_vec, cluster_rows = FALSE,
                        show_heatmap_legend = FALSE, height = unit(16.5, "cm"), width = unit(0.5, "cm"))
fig_upper <- hm_upper + inst_legend_upper +
  # add in names bc adding 2nd heatmap deletes them
  rowAnnotation(rn = anno_text(new_names_heatmap[colnames(upper_tert_reordered)]))
draw(fig_upper)

clusters_low <- hclust(d = dist(x = split.corr[[2]], method = "euclidean"), method = "complete")
lower_tert_reordered <- data.frame(split.corr[[2]]) %>%
  select(clusters_low[["order"]]) %>%
  t() %>%
  data.frame() %>%
  select(clusters_low[["order"]]) %>%
  as.matrix()
( hm_lower <- Heatmap(lower_tert_reordered, name = "Lower tertile z matrix", 
        col = col2, show_row_dend = FALSE, show_column_dend = FALSE,
        cluster_rows = FALSE, cluster_columns = FALSE,
        # human-readable names
        row_labels = new_names_heatmap[colnames(lower_tert_reordered)], column_labels = new_names_heatmap[colnames(lower_tert_reordered)],
        # square
        height = unit(16.5, "cm"), width = unit(16.5, "cm"),
        row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 7)) )
# now with instrument colors
inst_legend_lower <- Heatmap(as.matrix(colnames(lower_tert_reordered)), name = "Instrument Name",
                             col = inst_colors_vec, cluster_rows = FALSE,
                             show_heatmap_legend = FALSE, height = unit(16.5, "cm"), width = unit(0.5, "cm"))
fig_lower <- hm_lower + inst_legend_lower +
  # add in names bc adding 2nd heatmap deletes them
  rowAnnotation(rn = anno_text(new_names_heatmap[colnames(lower_tert_reordered)]))
draw(fig_lower)

#####################################################################
### COMPARE OBSERVED DIFFERENCE MATRICES TO THE NULL DISTRIBUTION ###
#####################################################################

# edit null_dist funtion from previous scripts bc now we want to randomly select two thirds of the data
NULL_DIST_Z <- function(x){   # input x is a df in which rows are observations, columns are subscales
  # randomly generate vector of row numbers that is double the length of the tertile splits
  sample_size <- floor(0.66*nrow(x))
  selected <- sample(x = nrow(x), size = sample_size)
  # then take the top and bottom half of this random vector of indices
  inda <- head(selected, n = (0.5 * length(selected)))
  indb <- tail(selected, n = (0.5 * length(selected)))
  # split df based on these indices
  terta <- x[inda, ]
  tertb <- x[indb, ]
  # correlate each separately
  acorr <- cor(terta, use = "pairwise.complete.obs")
  bcorr <- cor(tertb, use = "pairwise.complete.obs")
  # verify row and column names are in the same order
  assertthat::assert_that(identical(colnames(acorr), colnames(bcorr)))
  assertthat::assert_that(identical(rownames(acorr), rownames(bcorr)))
  # z-transform
  az <- fisherz(acorr)
  bz <- fisherz(bcorr)
  # subtract one matrix from the other to get difference matrix
  diff <- az - bz
  return(diff)
}
# iterate the null dist function 10000 times over the "subscales" dataframe
n <- 10000
null <- map(seq_len(n), ~ NULL_DIST_Z(x = subscales[, 17:86])) 
# give each df a name other than a number, so that I can use the df names as column names when I collapse the list
names(null) <- paste("V", seq_len(n), sep = "")

# check row/col names of observed and null matrices are all the same
assertthat::assert_that(identical(colnames(lmh_z), colnames(null[["V9"]])))  # random df from list
assertthat::assert_that(identical(rownames(lmh_z), rownames(null[["V84"]]))) # random df from list

# pull out the upper triangle of each null difference matrix and put into a df
null.vec <- map_dfc(null, function(d){
  output <- d[upper.tri(d, diag = FALSE)]
  return(output)                         
}) ## output will be a df in which rows are pairs of subscales and columns are each iteration of the sampling
# assign more descriptive row names, that tracks which row is which pair of observations
kept <- which(upper.tri(null[[1]], diag = FALSE), arr.ind = TRUE)   # indices of row and column names that were kept with upper.tri
rownames(null.vec) <- paste(dimnames(null[[1]])[[2]][kept[,2]], dimnames(null[[1]])[[1]][kept[,1]], sep = "_x_")
# do the same to the observed values // vectorize
observed.vec <- data.frame(vals = lmh_z[upper.tri(lmh_z, diag = FALSE)])
rownames(observed.vec) <- paste(dimnames(lmh_z)[[2]][kept[ , 2]], dimnames(lmh_z)[[1]][kept[ , 1]], sep = "_x_")
observed.vec$pair <- rownames(observed.vec)


## FIND THE PERCENTILE AT WHICH EACH OBSERVED VALUE FALLS ##
percentiles <- data.frame(pair = matrix(rep(NA, nrow(null.vec))), percentile = matrix(rep(NA, nrow(null.vec))))
for (i in 1:nrow(null.vec)){
  percentiles$percentile[i] <- ecdf(null.vec[i, ])(observed.vec$vals[i])
  percentiles$pair[i] <- rownames(observed.vec)[i]
}
# NOTE: percentile values will be slightly different with each time the code is run, because the randomized row groupings
# change each time; however, the observed value will not change

# save the null distribution of delta-r values to use with other wisc subscales
write.csv(null.vec, file = "./modified_data/null_delta_r_distribution_33_ztransform.csv")
# save the list of null dist matrices to use in the WSBM code
save(null, file = "./modified_data/null_dist_list_33_ztransform.RData")


## ORGANIZE INTO A TABLE ##

# arrange in order of percentile at which the value falls
percentiles <- percentiles %>%
  # add in the p-value, in addition to the percentile at which it falls, and the FDR adjusted p-values
  mutate(p_value = ifelse(test = percentile > 0.5, yes = 1 - percentile, no = percentile)) %>%
  rowwise() %>%
  # add in the observed difference between the pairwise correlation from the tertile split
  mutate(observed_diff = observed.vec$vals[observed.vec$pair %in% pair]) %>%
  ungroup() 
# table showing the significant pairs of subscales
percentiles %>%
  arrange(percentile) %>%
  filter(p_value <= 0.025) %>% nrow()   # percentile reaches significance in 241 edges
which(percentiles$p_value <= (0.025/2415)) %>% length()  # 3 survive Bonferroni

# save
write.csv(percentiles, file = "./modified_data/deltar_percentiles_WISC_FSIQ_33_ztransform.csv")

# Supplemental Data 1
upper_vec <- data.frame(corr = split.corr[[1]][upper.tri(split.corr[[1]], diag = FALSE)], 
                        edge = paste(dimnames(split.corr[[1]])[[2]][kept[ , 2]], dimnames(split.corr[[1]])[[1]][kept[ , 1]], sep = "_x_"))
lower_vec <- data.frame(corr = split.corr[[2]][upper.tri(split.corr[[2]], diag = FALSE)], 
                        edge = paste(dimnames(split.corr[[2]])[[2]][kept[ , 2]], dimnames(split.corr[[2]])[[1]][kept[ , 1]], sep = "_x_"))
# confirm all in same order
assertthat::assert_that(identical(percentiles$pair, upper_vec$edge))
assertthat::assert_that(identical(percentiles$pair, lower_vec$edge))
# compile
sd1 <- percentiles  %>%
  mutate(upper_corr = upper_vec$corr, 
         lower_corr = lower_vec$corr) %>%
  select(pair, lower_corr, upper_corr, observed_diff, p_value)
# export
write.csv(sd1, file = "./modified_data/supplementaldata1.csv")

#########################################
### DIFF MAT WITH NS PAIRS GREYED OUT ###
#########################################

# ns matrix in same order as reordered matrix
lmh_ns <- lmh_reordered
# have to put the percentiles df in the order of the lower triangle of lmh
replace <- which(upper.tri(lmh_ns, diag = FALSE), arr.ind = TRUE)   # indices of row and column names that need to be replaced
order <- paste(dimnames(lmh_ns)[[2]][kept[,2]], dimnames(lmh_ns)[[1]][kept[,1]], sep = "_x_") # what the order will correspond to in the percentiles df
# the order of subscale names in some pairs is the reverse of the percentiles$pair names, so they are missed
missed <- order[which(order %in% percentiles$pair == FALSE)]
# swap order of subscales for missed edges
new_names <- as.data.frame(str_split_fixed(missed, "_x_", n = 2))
missed_fixed <- paste(new_names$V2, new_names$V1, sep = "_x_")
# indices of the edges that were missed
missed_idx <- which(order %in% missed == TRUE)
# replace missed edges with edges reordered
order[missed_idx] <- missed_fixed
# reorder percentiles_df so it matches the order of edges in the clustered/reordered difference matrices
percentiles_reordered <- percentiles[match(order, percentiles$pair), ]
# replace the observed diff with NA if the pairwise diff is NS
percentiles_reordered <- percentiles_reordered %>%
  mutate(observed_diff = ifelse(test = p_value > 0.025, yes = NA, no = observed_diff))
# replace the upper.tri with NAs if the pair is NS
lmh_ns[upper.tri(lmh_ns, diag = FALSE)] <- percentiles_reordered$observed_diff

# check they match 
assertthat::assert_that(all(lmh_reordered == lmh_ns, na.rm = TRUE))
assertthat::assert_that(identical(colnames(lmh_reordered), colnames(lmh_ns)))
assertthat::assert_that(identical(rownames(lmh_reordered), rownames(lmh_ns)))

# visualize
( hm_diff_ns <- Heatmap(lmh_ns, name = "Tertile difference matrix, !NS", col = fixedcolors,
        row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 7),
        # human-readable names
        row_labels = new_names_heatmap[colnames(lmh_ns)], column_labels = new_names_heatmap[colnames(lmh_ns)],
        # square
        height = unit(16.5, "cm"), width = unit(16.5, "cm"),
        # need to assign clustering functions using the complete results_mat without NAs for NS
        cluster_rows = FALSE, cluster_columns = FALSE) )
fig_diffmat_ns <- hm_diff_ns + inst_legend1 + 
  # add in names bc adding 2nd heatmap deletes them
  rowAnnotation(rn = anno_text(new_names_heatmap[colnames(lmh_ns)]))
draw(fig_diffmat_ns)

# save to compare
write.csv(lmh_ns, file = "./modified_data/diff_mat_WISC_FSIQ_ns_33.csv")
# avg of sig edges
mean(lmh_ns[upper.tri(lmh_ns, diag = FALSE)], na.rm = TRUE)
# range
range(lmh_ns[upper.tri(lmh_ns, diag = FALSE)], na.rm = TRUE)


#####################################################
### CORRELATION BETWEEN COHEN'S D AND AVG DELTA-Z ###
#####################################################

# avg delta-z for each subscale
avg_z <- apply(lmh_z, 1, mean, na.rm = TRUE)
# check same order
identical(names(avg_z), cohens_d$subscale)
# correlate
cor(avg_z, cohens_d$d)
cohens_d %>%
  # add average of delta-z for that row
  mutate(z = avg_z) %>%
  # plot
  ggplot(., aes(d, z)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Cohen's d", y = "Average ∆z") +
  geom_text(inherit.aes = FALSE, size = 7,
            aes(x = -0.75, y = 0.08,
                label = paste0("r = ", round(cor(avg_z, cohens_d$d), 3))))
  




##################################################################
### SUBSCALE X SUBSCALE SCATTERPLOTS IN HIGH AND LOW IQ GROUPS ###
##################################################################

subscales %>%
  mutate(group = case_when(
    WISC_FSIQ >= range(top33$WISC_FSIQ)[1] ~ "High IQ", 
    WISC_FSIQ <= range(bottom33$WISC_FSIQ)[2] ~ "Low IQ"                    
  )) %>%
  ggplot(., aes(SCARED_SR_Total, SDQ_EP)) +   # replace with any two subscales of interest
  geom_jitter(alpha = 0.3) +
  geom_smooth(method = "lm") +
  facet_wrap(~ group) +
  # add correlation coefficient 
  stat_fit_glance(method = "lm",
                  method.args = list(formula = y ~ x), 
                  label.y = "top", 
                  label.x = "right",
                  mapping = aes(label = sprintf('r~"="~%.3f',
                                                sqrt(stat(r.squared)))), 
                  parse = TRUE)



################################
### TABLE ONE - DEMOGRAPHICS ###
################################
# read in data on SES and race
ses <- read.csv("./iqanalyses_rawdata/9994_Barratt_20210310.csv") %>%
  # remove first row
  magrittr::extract(-1, ) %>%
  mutate(Barratt_Total = as.numeric(Barratt_Total))
demog <- read.csv("./iqanalyses_rawdata/9994_PreInt_Demos_Fam_20210310.csv") %>%
  magrittr::extract(-1, )
# add race and SES to subscales df
t1data <- left_join(subscales, ses[, c(1,19)], by = "Anonymized.ID") %>%
  left_join(., demog[, c(1, 12:14)], by = "Anonymized.ID") %>%
  # some rows are duplicated, so remove duplicate rows
  distinct() %>%
  # sex as factor
  mutate(Sex = as.factor(Sex)) %>%
  # replace numbers with char for race
  mutate(Child_Race = case_when(
    Child_Race == "0" ~ "White", 
    Child_Race == "1" ~ "Black", 
    Child_Race == "2" ~ "Hispanic", 
    Child_Race == "3" ~ "Asian", 
    Child_Race == "4" ~ "Indian", 
    Child_Race == "5" ~ "Native American Indian", 
    Child_Race == "6" ~ "American Indian/Alaskan Native", 
    Child_Race == "7" ~ "Native Hawaiian/Other Pacific Islander", 
    Child_Race == "8" ~ "Two or more races", 
    Child_Race == "9" ~ "Other race", 
    Child_Race == "10" ~ "Unknown"
  )) %>%
  # add var for IQ group
  filter(!is.na(WISC_FSIQ)) %>%
  mutate(group = case_when(
    WISC_FSIQ >= range(top33$WISC_FSIQ)[1] ~ "highIQ", 
    WISC_FSIQ <= range(bottom33$WISC_FSIQ)[2] ~ "lowIQ", 
    WISC_FSIQ < range(top33$WISC_FSIQ)[1] & WISC_FSIQ > range(bottom33$WISC_FSIQ)[2] ~ "middle"
  )) 

# if no race data, mark as "Not Reported"
t1data$Child_Race[is.na(t1data$Child_Race)] <- "Not Reported"
# variables I want in the table
t1vars <- colnames(t1data[c(11, 12, 92, 95, 93)])
# make Table 1
table1 <- CreateTableOne(data = t1data, vars = t1vars, 
                         strata = "group",  # give separate descriptive stats for high and low IQ group
                         addOverall = TRUE) # keep column with whole-sample data
print(table1) # will want to input this into R Markdown file to get prettier output

# export --> format in Excel
t1 <- print(table1) %>%
  data.frame() %>%
  rownames_to_column(var = "X") 
write.csv(t1, file = "./figures/table1data_33.csv")

# chi-sq to check prop of females in each tertile
chisq_mat <- as.matrix(cbind(table(top33$Sex), table(bottom33$Sex)))
colnames(chisq_mat) <- c("highIQ", "lowIQ")

chisq.test(x = chisq_mat) # x2 = 0.415, df = 1, p-val = 0.519 non sig




##################################################
### SUPP TABLE TWO - DESCRIPTIVE STATS FOR ALL VARS ###
##################################################
# add var for IQ group
ts2data <- subscales %>%
  filter(!is.na(WISC_FSIQ)) %>%
  mutate(group = case_when(
    WISC_FSIQ >= range(top33$WISC_FSIQ)[1] ~ "highIQ", 
    WISC_FSIQ <= range(bottom33$WISC_FSIQ)[2] ~ "lowIQ", 
    WISC_FSIQ < range(top33$WISC_FSIQ)[1] & WISC_FSIQ > range(bottom33$WISC_FSIQ)[2] ~ "middle"
  )) 
# check assigned groups correctly
assertthat::assert_that(all(ts2data$Anonymized.ID[ts2data$group == "lowIQ"] %in% bottom33$Anonymized.ID))
assertthat::assert_that(all(ts2data$Anonymized.ID[ts2data$group == "highIQ"] %in% top33$Anonymized.ID))
# variables I want in the table
ts2vars <- colnames(subscales[17:92])
# make Table 1
tables2 <- CreateTableOne(data = ts2data, vars = ts2vars, 
                         strata = "group", addOverall = TRUE) # give separate descriptive stats for high and low IQ group, and whole sample
print(tables2) 

# export --> edit in Excel
ts2 <- print(tables2) %>%
  data.frame() %>%
  rownames_to_column(var = "subscale") %>%
  # don't need to say mean/SD in every row
  mutate(subscale = str_remove_all(string = subscale, pattern = " \\(mean \\(SD\\)\\)")) %>%
  # want to replace p-values with the p-values generated from the null distribution
  mutate(cohens_d = NA,
         diff_means = NA,
         p_from_null = NA) %>%
  select(!c(p, test))
assertthat::assert_that(identical(ts2$subscale[2:71], observed_means$subscale))
# enter in info on where obs value falls, and diff in means/cohen's d
ts2$cohens_d[2:71] <- round(observed_means$cohens_d, 3)
ts2$diff_means[2:71] <- round(observed_means$difference, 3)
ts2$p_from_null[2:71] <- round(observed_means$p_value, 3)
# save
write.csv(ts2, file = "./figures/tablesupp2.csv") # Table S2 in paper



#########################################
### COMPARE SUBSAMPLE TO WHOLE SAMPLE ###   
#########################################

# behavioral data for all participants
alldat <- read.csv("./modified_data/ninth_release_behavioral_data2.csv")
# add in ses and race
alldat <- left_join(alldat, ses[, c(1,19)], by = "Anonymized.ID") %>%
  left_join(., demog[, c(1, 12:14)], by = "Anonymized.ID") %>%
  # some rows are duplicated, so remove duplicate rows
  distinct() %>%
  # sex as factor
  mutate(Sex = as.factor(Sex)) %>%
  # replace numbers with char for race
  mutate(Child_Race = case_when(
    Child_Race == "0" ~ "White", 
    Child_Race == "1" ~ "Black", 
    Child_Race == "2" ~ "Hispanic", 
    Child_Race == "3" ~ "Asian", 
    Child_Race == "4" ~ "Indian", 
    Child_Race == "5" ~ "Native American Indian", 
    Child_Race == "6" ~ "American Indian/Alaskan Native", 
    Child_Race == "7" ~ "Native Hawaiian/Other Pacific Islander", 
    Child_Race == "8" ~ "Two or more races", 
    Child_Race == "9" ~ "Other race", 
    Child_Race == "10" ~ "Unknown"
  ))

# take only the participants NOT in our subsample
excluded <- alldat %>% 
  filter(Anonymized.ID %in% subscales$Anonymized.ID == FALSE)

# compare to subsample
t.test(t1data$Age, excluded$Age)                        # age SIG -- our subsample is younger
chisq_mat2 <- as.matrix(cbind(table(t1data$Sex), table(excluded$Sex)))
colnames(chisq_mat2) <- c("subsample", "total")
chisq.test(x = chisq_mat2)                              # sex NOT SIG
t.test(t1data$Barratt_Total, excluded$Barratt_Total)    # ses SIG -- higher ses in our subsample
subsample_races <- as.vector(table(t1data$Child_Race))[c(2:5, 8:11)] # vector of number of children of each race in subsample
excluded_races <- as.vector(table(excluded$Child_Race)) / sum(as.vector(table(excluded$Child_Race))) # prop of children of each race in excluded sample
chisq.test(x = subsample_races, p = excluded_races)     # race SIG -- proportions differ





################################
################################
##### MISSINGNESS ANALYSES #####
################################
################################

# upper tri corr between all data and a subset with no missing data
z_mat_full <- subscales %>% 
  select(CBCL_AD_T : ICU_SR_Tot) %>%
  cor(., use = "pairwise.complete.obs") %>%
  fisherz()
ut_full <- z_mat_full[upper.tri(z_mat_full, diag = FALSE)]
z_mat_no_miss <- subscales %>%
  na.omit() %>% 
  select(CBCL_AD_T : ICU_SR_Tot) %>%
  cor(., use = "pairwise.complete.obs") %>%
  fisherz()
ut_no_miss <- z_mat_no_miss[upper.tri(z_mat_no_miss, diag = FALSE)]
# correlate
cor(ut_full, ut_no_miss) # r = 0.989

### CORRELATE UT OF # OF PAIRWISE OBS WITH UT OF DELTA-Z MATRIX
# matrix in which each cell records the number of pairwise observations
number_pairwise_obs <- psych::pairwiseCount(x = subscales[, 17:86], diagonal = FALSE)
assertthat::assert_that(identical(colnames(lmh), colnames(number_pairwise_obs)))
assertthat::assert_that(identical(rownames(lmh), rownames(number_pairwise_obs)))
# pull out upper triangles
ut_number <- number_pairwise_obs[upper.tri(number_pairwise_obs, diag = FALSE)]
ut_delta <- lmh_z[upper.tri(lmh_z, diag = FALSE)]
#correlate
cor(ut_number, ut_delta) # r=0.031, very low
plot(ut_number, ut_delta)



######################################
### TESTING FOR NON-LINEAR EFFECTS ###
######################################
# model each subscale based on IQ + IQ^2
non_lin <- list()
for(i in colnames(subscales[17:86])){
  non_lin[[i]] <- summary(lm(get(i) ~ WISC_FSIQ + I(WISC_FSIQ^2), subscales)) %>% 
    magrittr::extract2("coefficients") %>%
    data.frame()
}
# p-value of IQ^2 term is in the 3rd row, 4th column, so extract that value
non_lin_df <- map_dfc(.x = non_lin, .f = function(x){
  output <- x[3,4]
  return(output)
}) %>%
  t() %>%
  data.frame() %>%
  magrittr::set_colnames("iq2_term_pval")
# find how many are sig
non_lin_df %>%
  filter(iq2_term_pval <= (0.05/70))   # 6/70 subscales significant after Bonferroni correction

# df with beta values to use in polynomial fit equation in plots
sig_names <- non_lin_df %>%
  filter(iq2_term_pval <= (0.05/70)) %>%
  rownames()
betas <- data.frame()
for(n in sig_names){
  betas[n, 1] <- non_lin[[n]][["Estimate"]][1]
  betas[n, 2] <- non_lin[[n]][["Estimate"]][2]
  betas[n, 3] <- non_lin[[n]][["Estimate"]][3]
}
colnames(betas) <- c("y_int", "iq", "iq_sq")

# initialize list for plots
polynomial_plts <- list()
# function to plot polynomial fit line -- SRS_COG_T
f_srs_cog <- function(x) betas["SRS_COG_T",1] + betas["SRS_COG_T",2]*(x) +betas["SRS_COG_T",3]*(x^2)
# plot
polynomial_plts[[1]] <- ggplot(subscales, aes(WISC_FSIQ, SRS_COG_T)) +
  geom_point(size = 0.5) +
  geom_function(fun= f_srs_cog, color = "blue") +
  geom_smooth(method = "gam", color = "green", se = FALSE) +
  geom_smooth(method = "lm", color = "red", se = FALSE)
# function to plot polynomial fit line -- SRS_COM_T
f_srs_com <- function(x) betas["SRS_COM_T",1] + betas["SRS_COM_T",2]*(x) +betas["SRS_COM_T",3]*(x^2)
# plot
polynomial_plts[[2]] <- ggplot(subscales, aes(WISC_FSIQ, SRS_COM_T)) +
  geom_point(size = 0.5) +
  geom_function(fun= f_srs_com, color = "blue") +
  geom_smooth(method = "gam", color = "green", se = FALSE) +
  geom_smooth(method = "lm", color = "red", se = FALSE)
# function to plot polynomial fit line -- SRS_SCI_T
f_srs_sci <- function(x) betas["SRS_SCI_T",1] + betas["SRS_SCI_T",2]*(x) +betas["SRS_SCI_T",3]*(x^2)
# plot
polynomial_plts[[3]] <- ggplot(subscales, aes(WISC_FSIQ, SRS_SCI_T)) +
  geom_point(size = 0.5) +
  geom_function(fun= f_srs_sci, color = "blue") +
  geom_smooth(method = "gam", color = "green", se = FALSE) +
  geom_smooth(method = "lm", color = "red", se = FALSE)
# function to plot polynomial fit line -- SRS_Total_T
f_srs_tot <- function(x) betas["SRS_Total_T",1] + betas["SRS_Total_T",2]*(x) +betas["SRS_Total_T",3]*(x^2)
# plot
polynomial_plts[[4]] <- ggplot(subscales, aes(WISC_FSIQ, SRS_Total_T)) +
  geom_point(size = 0.5) +
  geom_function(fun= f_srs_tot, color = "blue") +
  geom_smooth(method = "gam", color = "green", se = FALSE) +
  geom_smooth(method = "lm", color = "red", se = FALSE)
# function to plot polynomial fit line -- RBS_RESTR_INTERESTS
f_rbs_ri <- function(x) betas["RBS_RI",1] + betas["RBS_RI",2]*(x) +betas["RBS_RI",3]*(x^2)
# plot
polynomial_plts[[5]] <- ggplot(subscales, aes(WISC_FSIQ, RBS_RI)) +
  geom_point(size = 0.5) +
  geom_function(fun= f_rbs_ri, color = "blue") +
  geom_smooth(method = "gam", color = "green", se = FALSE) +
  geom_smooth(method = "lm", color = "red", se = FALSE)
# function to plot polynomial fit line -- RBS_Total
f_rbs_tot <- function(x) betas["RBS_Tot",1] + betas["RBS_Tot",2]*(x) +betas["RBS_Tot",3]*(x^2)
# plot
polynomial_plts[[6]] <- ggplot(subscales, aes(WISC_FSIQ, RBS_Tot)) +
  geom_point(size = 0.5) +
  geom_function(fun= f_rbs_tot, color = "blue") +
  geom_smooth(method = "gam", color = "green", se = FALSE) +
  geom_smooth(method = "lm", color = "red", se = FALSE)

# arrange all plots together
ggarrange(plotlist = polynomial_plts)


### ORTHOGONAL POLYNOMIAL CODING ###
subscales_ord <- subscales %>%
  mutate(iq_ord = as.factor(case_when(
    WISC_FSIQ <= 91 ~ "low", 
    WISC_FSIQ < 106 & WISC_FSIQ > 91 ~ "mid", 
    WISC_FSIQ >= 106 ~ "high"
  ))) 
# confirm 937 low and 936 high
length(which(subscales_ord$iq_ord == "low"))
length(which(subscales_ord$iq_ord == "high"))
# set low < mid < high in order
subscales_ord$iq_ord <- ordered(subscales_ord$iq_ord, levels = c("low", "mid", "high"))

contrasts(subscales_ord$iq_ord) = contr.poly(3)
summary(lm(SRS_COG_T ~ iq_ord, subscales_ord)) # sig linear
summary(lm(SRS_COM_T ~ iq_ord, subscales_ord)) # sig linear
summary(lm(SRS_SCI_T ~ iq_ord, subscales_ord)) # sig linear
summary(lm(SRS_Total_T ~ iq_ord, subscales_ord)) # sig linear
summary(lm(RBS_RI ~ iq_ord, subscales_ord)) # sig linear
summary(lm(RBS_Tot ~ iq_ord, subscales_ord)) # sig linear
# ALL LINEAR TERMS ARE SIGNIFICANT, NO SIGNIFICANT QUADRATIC TERMS


#################################################
### SUPPLEMENTAL D: HCLUST/CUTREE ON TERTILES ###
#################################################

### function for generating df for the step-wise height x k plot ###
HEIGHT_X_K <- function(dendrogram, start_height, end_height, interval){
  # create df to enter height and number of groupings (k) that result from that cut
  dend_groupings <- data.frame(height = seq(from = start_height, to = end_height, by = interval), 
                               k = rep(NA))
  # cut tree at different heights, in intervals of 0.01 
  for (i in seq(from = start_height, to = end_height, by = interval)){
    # cut at each interval
    dend_cut <- cutree(dendrogram, h = i) 
    # calculate number of groups and enter into df
    dend_groupings[dend_groupings$height == i, 2] <- length(unique(dend_cut))
  }
  return(dend_groupings)
}
# build dendrograms
dend_top33 <- hclust(dist(x = cor(top33[, 17:86], use = "pairwise.complete.obs"), method = "euclidean"), method = "complete")
dend_bottom33 <- hclust(dist(x = cor(bottom33[, 17:86], use = "pairwise.complete.obs"), method = "euclidean"), method = "complete")
# plot
plot(dend_top33)
plot(dend_bottom33)

# generate dfs of height x k clusters
dend_groupings_top33 <-  HEIGHT_X_K(dend_top33, start_height = 0, end_height = 4.5, interval = 0.01)
dend_groupings_bottom33 <-  HEIGHT_X_K(dend_bottom33, start_height = 0, end_height = 4.5, interval = 0.01)
# merge
dend_groupings_tertiles <- data.frame(height = dend_groupings_top33$height, 
                                      upper_tertile_k = dend_groupings_top33$k, 
                                      lower_tertile_k = dend_groupings_bottom33$k) %>%
  # long form
  pivot_longer(cols = !height, names_to = "tertile", values_to = "k")
# plot
ggplot(dend_groupings_tertiles, aes(height, k)) +
  geom_step(aes(color = tertile)) +
  # can truncate the y axis so that it zooms in on the most important divisions
  coord_cartesian(ylim = c(0, 30)) +
  # have the steps proceed in the other direction
  scale_x_reverse() +
  # clean up labels
  labs(y = "Number of clusters, k", x = "Cut height", color = "Tertile") +
  scale_color_discrete(labels = c("Lower tertile", "Upper tertile"))


############################################################
### RESAMPLE TO GENERATE SMOOTH LINE/CONFIDENCE INTERVAL ###
############################################################

BOOTSTRAP_FXN <- function(upper, lower){ # upper is upper tertile, lower is lower tertile, n is # of iterations
  # resample each tertile with replacement
  a <- slice_sample(lower, prop = 1, replace = TRUE)
  b <- slice_sample(upper, prop = 1, replace = TRUE)
  # build dendrograms
  dend_a <- hclust(dist(x = cor(a[, 17:86], use = "pairwise.complete.obs"), method = "euclidean"), method = "complete")
  dend_b <- hclust(dist(x = cor(b[, 17:86], use = "pairwise.complete.obs"), method = "euclidean"), method = "complete")
  # generate dfs of height x k clusters
  dend_groupings_a <-  HEIGHT_X_K(dend_a, start_height = 0, end_height = 4.5, interval = 0.01)
  dend_groupings_b <-  HEIGHT_X_K(dend_b, start_height = 0, end_height = 4.5, interval = 0.01)
  # merge into one df
  dend_groupings_ab <- data.frame(height = dend_groupings_a$height, 
                                  upper_tertile_k = dend_groupings_b$k, 
                                  lower_tertile_k = dend_groupings_a$k)
  return(dend_groupings_ab)
}
# run this function 1000 times to generate CI
ab_results <- map(seq_len(1000), ~ BOOTSTRAP_FXN(upper = top33, lower = bottom33)) 
names(ab_results) <- paste("V", seq_len(1000), sep = "")
for (i in 1:1000){
  ab_results[[i]]$iteration <- paste("V", i, sep = "")
}
# collapse into one df
ab_results_df <- Reduce(rbind, ab_results)
# put in long form
ab_results_df_long <- pivot_longer(ab_results_df, cols = !c(iteration, height), names_to = "tertile", values_to = "k")
# plot
ggplot(ab_results_df_long, aes(height, k)) +
  geom_step(aes(color = tertile), alpha = 0.1) +
  # can truncate the y axis so that it zooms in on the most important divisions
  coord_cartesian(ylim = c(0, 30)) +
  # have the steps proceed in the other direction
  scale_x_reverse()

# average from resampling
ab_avg <- ab_results_df %>%
  group_by(height) %>%
  summarise(avg_lower = mean(lower_tertile_k), 
            avg_upper = mean(upper_tertile_k)) %>%
  ungroup()
# plot average lines
ab_avg %>%
  pivot_longer(cols = !height, names_to = "tertile", values_to = "avg_k") %>%
  ggplot(., aes(height, avg_k)) +
  geom_line(aes(color = tertile)) +
  coord_cartesian(ylim = c(0, 30)) +
  # have the steps proceed in the other direction
  scale_x_reverse() +
  # clean up labels
  labs(y = "Average number of clusters, k", x = "Cut height", color = "Tertile") +
  scale_color_discrete(labels = c("Lower tertile", "Upper tertile"))

# 95th and 5th centile lines
ci_step <- ab_results_df_long %>%
  group_by(height) %>% 
  mutate(perc95_upper = quantile(k[tertile == "upper_tertile_k"], probs = 0.95), 
         perc5_upper = quantile(k[tertile == "upper_tertile_k"], probs = 0.05), 
         perc95_lower = quantile(k[tertile == "lower_tertile_k"], probs = 0.95), 
         perc5_lower = quantile(k[tertile == "lower_tertile_k"], probs = 0.05)) %>%
  select(height, starts_with("perc")) %>%
  # only need the unique rows
  magrittr::extract(seq(1, 902, 2), )

# add in observed values
dend_groupings_tertiles2 <- dend_groupings_tertiles %>%
  pivot_wider(names_from = "tertile", values_from = "k")
ci_step %>%
  rowwise() %>%
  mutate(upper_observed = dend_groupings_tertiles2$upper_tertile_k[dend_groupings_tertiles2$height %in% height], 
         lower_observed = dend_groupings_tertiles2$lower_tertile_k[dend_groupings_tertiles2$height %in% height]) %>%
  ungroup() %>%
  pivot_longer(cols = !c(height), names_to = "percentile", values_to = "k") %>%
  # plot
  ggplot(., aes(height, k)) +
  geom_step(aes(color = percentile)) +
  # can truncate the y axis so that it zooms in on the most important divisions
  coord_cartesian(ylim = c(0, 30)) +
  # have the steps proceed in the other direction
  scale_x_reverse() +
  # set color values
  scale_color_manual(values = c("red", "darkorange", "powderblue", "darkorange", "powderblue", "blue"))



#####################################
##### SPLINE REGRESSION AGE/SEX #####
#####################################

# load in df with raw scores only, none t-scored
raw_subscales <- read.csv("./modified_data/ninth_release_rawscores2.csv")
# sex as ordinal var
raw_subscales <- raw_subscales %>%
  mutate(sex.ord = as.ordered(Sex))
# apply gam to each behavioral subscale raw score, to regress out age and sex
model1 <- apply(raw_subscales[17:86], 2, FUN = function(x){
  gam(data = raw_subscales, x ~ s(raw_subscales$Age) + raw_subscales$sex.ord + s(raw_subscales$Age, by = raw_subscales$sex.ord))
})

# extract residuals
model1_residuals <- map(model1, function(x) residuals(x)) 
# enter residuals into df where raw scores were previously--need to be aware that gam() removed NA values
agesex_residuals <- data.frame(matrix(NA, nrow = 2753, ncol = 70))
colnames(agesex_residuals) <- colnames(raw_subscales[17:86])
for (n in colnames(agesex_residuals)){
  temp <- as.vector(raw_subscales[, n])
  resid <- model1_residuals[[n]]
  temp[is.na(temp) == FALSE] <- resid
  agesex_residuals[, n] <- temp
}
# add in rest of the info
agesex_residuals[, 71:92] <- raw_subscales[, c(1:16, 87:92)]

# split into IQ tertiles
top33_bas <- agesex_residuals %>%
  filter(!is.na(WISC_FSIQ)) %>%    # remove obs with no WISC IQ data
  slice_max(., WISC_FSIQ, prop = 0.33)
bottom33_bas <- agesex_residuals %>%
  filter(!is.na(WISC_FSIQ)) %>%    # remove obs with no WISC IQ data
  slice_min(., WISC_FSIQ, prop = 0.33)
# collapse into a list
split_bas <- list(high = top33_bas, low = bottom33_bas)
# check that same participants are in top and bottom 30 as in raw dataset
assertthat::assert_that(identical(top33_bas$Anonymized.ID, top33$Anonymized.ID))
assertthat::assert_that(identical(bottom33_bas$Anonymized.ID, bottom33$Anonymized.ID))



########################################
### COHEN'S D WITH AGE/SEX RESIDUALS ###
########################################

# calculate means
means_bas <- map_dfc(split_bas, function(d){
  d %>%
    # select subscales
    select(CBCL_AD : ICU_SR_Tot) %>%
    # compute mean
    apply(., 2, mean, na.rm = TRUE) %>%
    data.frame(.) %>%
    magrittr::set_rownames(colnames(d[1:70]))
}) 
# clean up
observed_means_bas <- means_bas %>%
  set_colnames(c("high", "low")) %>%
  rownames_to_column("subscale") %>%
  # add a var for the difference between mean score in low vs high IQ groups
  mutate(difference = high - low) 

# calculate standard deviation of each subscale in the high and low IQ groups
stdev_bas <- map_dfc(split_bas, function(d){
  d %>%
    # select subscales
    select(CBCL_AD : ICU_SR_Tot) %>%
    # compute standard deviation
    apply(., 2, sd, na.rm = TRUE) %>%
    # make df
    data.frame(.) %>%
    # keep subscales as row names
    magrittr::set_rownames(colnames(d[1:70]))
}) %>%
  # edits to colnames and save subscales as a column, not rownames
  set_colnames(c("high", "low")) %>%
  rownames_to_column("subscale") %>%
  # compute the pooled standard deviation, which will be used to caluclate cohen's d
  mutate(pooled_sd = sqrt(((high ^ 2) + (low ^ 2)) / 2) )

# divide the difference in means by the pooled standard deviation to get cohen's d
# first check that rows are in the same order
assertthat::assert_that(identical(observed_means_bas$subscale, stdev_bas$subscale))
cohens_d_bas <- data.frame(subscale = observed_means_bas$subscale, d = observed_means_bas$difference / stdev_bas$pooled_sd) 
# present in order of highest to lowest effect size
cohens_d_bas %>%
  arrange(desc(d)) 

# compare to cohen's d using mix of raw and normed scores
cor(cohens_d$d, cohens_d_bas$d) # r= 0.988; very high
# plot with identity line
cohens_d %>%
  select(subscale, d) %>%
  dplyr::rename(main = d) %>%
  mutate(bas = cohens_d_bas$d[1:70]) %>%
  ggplot(., aes(x = main, y = bas)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Cohen's d values in main analyses", y = "Cohen's d values controlling for age and sex")


##################################################
### CORR MATRIX (Z-TRANSFORMED) WITH RESIDUALS ###
##################################################
corr_bas_z <- cor(agesex_residuals[, 1:70], use = "pairwise.complete.obs") %>%
  fisherz()
corr_reg_z <- cor(subscales[17:86], use = "pairwise.complete.obs") %>%
  fisherz()
# compare-- checked rows and cols in same order
cor(corr_bas_z[upper.tri(corr_bas_z, dia = FALSE)], corr_reg_z[upper.tri(corr_reg_z, dia = FALSE)]) # r=0.993

##################################
### DIFF MATRIX WITH RESIDUALS ###
##################################
# calculate correlation matrix for each tertile, then z-transform
split.corr_bas <- map(split_bas, function(d){
  d %>% 
    select(CBCL_AD : ICU_SR_Tot) %>%
    cor(., use = "pairwise.complete.obs") %>%
    fisherz()
})

# subtract LOW from HIGH to make difference matrix
assertthat::assert_that(identical(colnames(split.corr_bas[["low"]]),  colnames(split.corr_bas[["high"]]))) # check row and col names are identical
assertthat::assert_that(identical(rownames(split.corr_bas[["low"]]),  rownames(split.corr_bas[["high"]])))
diffmat_bas <- split.corr_bas[["high"]] - split.corr_bas[["low"]] 

# compare
cor(diffmat_bas[upper.tri(diffmat_bas, dia = FALSE)], lmh_z[upper.tri(lmh_z, dia = FALSE)]) # r=0.945
#save
write.csv(diffmat_bas, file = "./modified_data/diff_mat_age_sex_regressed_z.csv")


### compare to null distribution
n=10000
null_bas <- map(seq_len(n), ~ NULL_DIST_Z(x = agesex_residuals[, 1:70])) 
# give each df a name other than a number, so that I can use the df names as column names when I collapse the list
names(null_bas) <- paste("V", seq_len(n), sep = "")

# check row/col names of observed and null matrices are all the same
assertthat::assert_that(identical(colnames(diffmat_bas), colnames(null_bas[["V9"]])))  # random df from list
assertthat::assert_that(identical(rownames(diffmat_bas), rownames(null_bas[["V84"]]))) # random df from list

# pull out the upper triangle of each null difference matrix and put into a df
null.vec_bas <- map_dfc(null_bas, function(d){
  output <- d[upper.tri(d, diag = FALSE)]
  return(output)                         
}) ## output will be a df in which rows are pairs of subscales and columns are each iteration of the sampling
# assign more descriptive row names, that tracks which row is which pair of observations
kept_bas <- which(upper.tri(null_bas[[1]], diag = FALSE), arr.ind = TRUE)   # indices of row and column names that were kept with upper.tri
rownames(null.vec_bas) <- paste(dimnames(null_bas[[1]])[[2]][kept_bas[,2]], dimnames(null_bas[[1]])[[1]][kept_bas[,1]], sep = "_x_")
# do the same to the observed values // vectorize
observed.vec_bas <- data.frame(vals = diffmat_bas[upper.tri(diffmat_bas, diag = FALSE)])
rownames(observed.vec_bas) <- paste(dimnames(diffmat_bas)[[2]][kept_bas[ , 2]], 
                                    dimnames(diffmat_bas)[[1]][kept_bas[ , 1]], sep = "_x_")
observed.vec_bas$pair <- rownames(observed.vec_bas)

## FIND THE PERCENTILE AT WHICH EACH OBSERVED VALUE FALLS ##
percentiles_bas <- data.frame(pair = matrix(rep(NA, nrow(null.vec_bas))), percentile = matrix(rep(NA, nrow(null.vec_bas))))
for (i in 1:nrow(null.vec_bas)){
  percentiles_bas$percentile[i] <- ecdf(null.vec_bas[i, ])(observed.vec_bas$vals[i])
  percentiles_bas$pair[i] <- rownames(observed.vec_bas)[i]
}
# arrange in order of percentile at which the value falls
percentiles_bas <- percentiles_bas %>%
  # add in the p-value, in addition to the percentile at which it falls, and the FDR adjusted p-values
  mutate(p_value = ifelse(test = percentile > 0.5, yes = 1 - percentile, no = percentile)) %>%
  rowwise() %>%
  # add in the observed difference between the pairwise correlation from the tertile split
  mutate(observed_diff = observed.vec_bas$vals[observed.vec_bas$pair %in% pair]) %>%
  ungroup() 
# table showing the significant pairs of subscales
percentiles_bas %>%
  arrange(percentile) %>%
  filter(p_value <= 0.025) %>% nrow()   # percentile reaches significance
which(percentiles_bas$p_value <= (0.025/2415)) %>% length()  # 4 survive Bonferroni

# save
write.csv(percentiles_bas, file = "./modified_data/deltar_percentiles_agesex_resid.csv")
# save iq resid null dist
write.csv(null.vec_bas, file = "./modified_data/null_delta_r_distribution_agesex_resid.csv")







###############################################################
### CONTROLLING FOR EFFECTS OF SES AND RACE/ETHNICITY ON IQ ###
###############################################################

# build linear model that residualized WISC_FSIQ using RACE and ses
linmod <- lm(WISC_FSIQ ~ as.numeric(Barratt_Total) + as.factor(Child_Race), data = t1data)
# now, make a df of all the subscale scores that only include the observations that were included in the lm (some were left out bc of incomplete data)
data_riq <- t1data[rownames(linmod[["model"]]), ] %>%
  # add the residualized IQ score column
  mutate(fsiq_resid = linmod[["residuals"]])
# check entered data in correct order; if not true stop the code
assertthat::assert_that(all(linmod[["model"]]["WISC_FSIQ"] == data_riq$WISC_FSIQ))

# dfs of top and bottom 33% of IQ distribution
top33_resid_iq <- data_riq %>%
  slice_max(., fsiq_resid, prop = 0.33)
bottom33_resid_iq <- data_riq %>%
  slice_min(., fsiq_resid, prop = 0.33)
# collapse into a list
split_resid_iq <- list(high = top33_resid_iq, low = bottom33_resid_iq)

###############################################
### IQ X SUBSCALE CORRELATIONS RESID FOR IQ ###
###############################################
# calculate iq x subscale correlations
correlations_riq <- cor(t1data[17:92], use = "pairwise.complete.obs")
cor_with_iq_resid <- data.frame(corr = correlations_riq[, 76]) %>%
  rownames_to_column(var = "subscale")
# correlate with uncorrected data
identical(cor_with_iq$subscale, cor_with_iq_resid$subscale)
cor(cor_with_iq$corr, cor_with_iq_resid$corr) # r = 1




####################################
### CORRELATING COHENS D VECTORS ###
####################################
means_riq <- map_dfc(split_resid_iq, function(d){
  d %>%
    # select subscales
    select(CBCL_AD_T : WISC_FSIQ) %>%
    # compute mean
    apply(., 2, mean, na.rm = TRUE) %>%
    data.frame(.) %>%
    magrittr::set_rownames(colnames(d[17:92]))
})
# clean up
observed_means_riq <- means_riq %>%
  set_colnames(c("high", "low")) %>%
  rownames_to_column("subscale") %>%
  # add a var for the difference between mean score in low vs high IQ groups
  mutate(difference = high - low) 

# calculate standard deviation of each subscale in the high and low IQ groups
stdev_riq <- map_dfc(split_resid_iq, function(d){
  d %>%
    # select subscales
    select(CBCL_AD_T : WISC_FSIQ) %>%
    # compute standard deviation
    apply(., 2, sd, na.rm = TRUE) %>%
    # make df
    data.frame(.) %>%
    # keep subscales as row names
    magrittr::set_rownames(colnames(d[17:92]))
}) %>%
  # edits to colnames and save subscales as a column, not rownames
  set_colnames(c("high", "low")) %>%
  rownames_to_column("subscale") %>%
  # compute the pooled standard deviation, which will be used to caluclate cohen's d
  mutate(pooled_sd = sqrt(((high ^ 2) + (low ^ 2)) / 2) )

# divide the difference in means by the pooled standard deviation to get cohen's d
# first check that rows are in the same order
assertthat::assert_that(identical(observed_means_riq$subscale, stdev_riq$subscale))
cohens_d_riq <- data.frame(subscale = observed_means_riq$subscale, d = observed_means_riq$difference / stdev_riq$pooled_sd) 
# present in order of highest to lowest effect size
cohens_d_riq %>%
  arrange(desc(d)) 

# compare to cohen's d splitting on non-corrected IQ scores
identical(cohens_d$subscale, cohens_d_riq$subscale[1:70])
cor(cohens_d$d, cohens_d_riq$d[1:70]) # r= 0.939; very high

# plot with identity line
cohens_d %>%
  select(subscale, d) %>%
  dplyr::rename(main = d) %>%
  mutate(riq = cohens_d_riq$d[1:70]) %>%
  ggplot(., aes(x = main, y = riq)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Cohen's d values in main analyses", y = "Cohen's d values controlling for SES and race")

####################################
### CORRELATING DELTA-Z MATRICES ###
####################################
# calculate correlations
split.corr_riq <- map(split_resid_iq, function(d){
  d %>% 
    select(CBCL_AD_T : ICU_SR_Tot) %>%
    cor(., use = "pairwise.complete.obs") %>%
    # then z-transform
    fisherz()
})
# subtract high from low to make difference matrix, in which positive difference means that pair of subscales are more highly
# correlated in high IQ than low IQ individuals
identical(colnames(split.corr_riq[["low"]]),  colnames(split.corr_riq[["high"]])) # check row and col names are identical
identical(rownames(split.corr_riq[["low"]]),  rownames(split.corr_riq[["high"]]))
diffmat_riq <- split.corr_riq[["high"]] - split.corr_riq[["low"]]               # subtract

# compare to non-normed iq split results
identical(colnames(lmh_z), colnames(diffmat_riq))
identical(rownames(lmh_z), rownames(diffmat_riq))
cor(diffmat_riq[upper.tri(diffmat_riq, diag = FALSE)], lmh_z[upper.tri(lmh_z, diag = FALSE)]) # r=0.873
#SAVE
write.csv(diffmat_riq, file = "./modified_data/diff_mat_iq_regressed_z.csv")


# compare to null distribution
null_riq <- map(seq_len(n), ~ NULL_DIST_Z(x = data_riq[, 17:86])) 
# give each df a name other than a number, so that I can use the df names as column names when I collapse the list
names(null_riq) <- paste("V", seq_len(n), sep = "")

# check row/col names of observed and null matrices are all the same
assertthat::assert_that(identical(colnames(diffmat_riq), colnames(null_riq[["V9"]])))  # random df from list
assertthat::assert_that(identical(rownames(diffmat_riq), rownames(null_riq[["V84"]]))) # random df from list

# pull out the upper triangle of each null difference matrix and put into a df
null.vec_riq <- map_dfc(null_riq, function(d){
  output <- d[upper.tri(d, diag = FALSE)]
  return(output)                         
}) ## output will be a df in which rows are pairs of subscales and columns are each iteration of the sampling
# assign more descriptive row names, that tracks which row is which pair of observations
kept_riq <- which(upper.tri(null_riq[[1]], diag = FALSE), arr.ind = TRUE)   # indices of row and column names that were kept with upper.tri
rownames(null.vec_riq) <- paste(dimnames(null_riq[[1]])[[2]][kept_riq[,2]], dimnames(null_riq[[1]])[[1]][kept_riq[,1]], sep = "_x_")
# do the same to the observed values // vectorize
observed.vec_riq <- data.frame(vals = diffmat_riq[upper.tri(diffmat_riq, diag = FALSE)])
rownames(observed.vec_riq) <- paste(dimnames(diffmat_riq)[[2]][kept_riq[ , 2]], 
                                    dimnames(diffmat_riq)[[1]][kept_riq[ , 1]], sep = "_x_")
observed.vec_riq$pair <- rownames(observed.vec_riq)

## FIND THE PERCENTILE AT WHICH EACH OBSERVED VALUE FALLS ##
percentiles_riq <- data.frame(pair = matrix(rep(NA, nrow(null.vec_riq))), percentile = matrix(rep(NA, nrow(null.vec_riq))))
for (i in 1:nrow(null.vec_riq)){
  percentiles_riq$percentile[i] <- ecdf(null.vec_riq[i, ])(observed.vec_riq$vals[i])
  percentiles_riq$pair[i] <- rownames(observed.vec_riq)[i]
}
# arrange in order of percentile at which the value falls
percentiles_riq <- percentiles_riq %>%
  # add in the p-value, in addition to the percentile at which it falls, and the FDR adjusted p-values
  mutate(p_value = ifelse(test = percentile > 0.5, yes = 1 - percentile, no = percentile), 
         adjusted_p_value = p.adjust(p_value, method = "fdr")) %>%
  rowwise() %>%
  # add in the observed difference between the pairwise correlation from the tertile split
  mutate(observed_diff = observed.vec_riq$vals[observed.vec_riq$pair %in% pair]) %>%
  ungroup() 
# table showing the significant pairs of subscales
percentiles_riq %>%
  arrange(percentile) %>%
  filter(p_value <= 0.025) %>% nrow()   # percentile reaches significance 186
which(percentiles_riq$p_value <= (0.025/2415)) %>% length()  # 2 survive Bonferroni
# save
write.csv(percentiles_riq, file = "./modified_data/deltar_percentiles_iq_resid.csv")
# save iq resid null dist
write.csv(null.vec_riq, file = "./modified_data/null_delta_r_distribution_iq_resid.csv")




#################################################
#################################################
### ANALYSES IN RESPONSE TO REVIEWER COMMENTS ###
#################################################
#################################################

# load in Dx data from clinician consensus
cc <- read.csv("raw data/9994_ConsensusDx_20210310.csv")[-1,]

# only include participants who are included in our subset with IQ data
cc_sub <- cc[cc$Anonymized.ID %in% subscales$Anonymized.ID, ] # n = 4605 so there must be some duplicates

# distinct() takes rows which are different from all others in at least one column (i.e., if 9/10 columns
# are the same but the 10th is different, the two rows are distinct) 
distinct_cc <- distinct(cc_sub)
# duplicate_distinct will show which distinct rows contain the same anonymized id, to try to determine
# what's different about them
duplicate_distinct <- distinct_cc[distinct_cc$Anonymized.ID %in% distinct_cc[which(duplicated(distinct_cc$Anonymized.ID) == TRUE), "Anonymized.ID"], ] %>%
  arrange(Anonymized.ID)

# check if the visit year and season is the same for each duplicate pair
duplicate_distinct %>% 
  group_by(Anonymized.ID, Year, Season) %>%
  filter(n() != 2) 
# one participant, A00093547, year/season are different from each other and from what's listed in subscales df 

# same idea, for Dx categories
test <- duplicate_distinct %>% 
  group_by(Anonymized.ID, DX_01_Cat) %>%
  filter(n() != 2) # 25 participants have different entries for Dx_01_Cat; most of them are differences in coding "disruptive" vs "disruptive, impulse control, and conduct"

test2 <- duplicate_distinct %>% 
  group_by(Anonymized.ID, DX_02_Cat) %>%
  filter(n() != 2) # 25 participants have different entries for Dx_01_Cat; most of them are differences in coding "disruptive" vs "disruptive, impulse control, and conduct"

# what sub-Dx are listed under NDD category?
cc %>%
  filter(DX_02_Cat == "Neurodevelopmental Disorders") %>%
  select(DX_02_Sub) %>%
  unique() # spec. learning disorder, ADHD, ASD, IntDis, communication disorder, motor disorder, other

# removing all participants who are duplicated, bc can't be sure which entry is accurate
cc_final <- distinct_cc[!(distinct_cc$Anonymized.ID %in% duplicate_distinct$Anonymized.ID), ]

# grouping by category (using _Cat for all except NDDs, for which using _Sub)
adhd <- cc_final %>%
  filter_all(any_vars(str_detect(., 'Attention-Deficit/Hyperactivity Disorder')))
asd <- cc_final %>%
  filter_all(any_vars(str_detect(., 'Autism Spectrum Disorder')))
anxiety <- cc_final %>%
  filter_all(any_vars(str_detect(., 'Anxiety Disorders')))
disruptive <- cc_final %>%
  filter_all(any_vars(str_detect(., 'Disruptive, Impulse Control and Conduct Disorders')))
depressive <-  cc_final %>%
  filter_all(any_vars(str_detect(., 'Depressive Disorders')))
learning <- cc_final %>%
  filter_all(any_vars(str_detect(., 'Specific Learning Disorder')))
communication <- cc_final %>%
  filter_all(any_vars(str_detect(., 'Communication Disorder')))
trauma <- cc_final %>%
  filter_all(any_vars(str_detect(., 'Trauma and Stressor Related Disorders')))
intellectual_dis <- cc_final %>%
  filter_all(any_vars(str_detect(., 'Intellectual Disability')))
motor_dis <- cc_final %>%
  filter_all(any_vars(str_detect(., 'Motor Disorder')))
ocd <- cc_final %>%
  filter_all(any_vars(str_detect(., "Obsessive Compulsive and Related Disorders")))
scz <- cc_final %>%
  filter_all(any_vars(str_detect(., "Schizophrenia Spectrum and other Psychotic Disorders")))
neurocog <- cc_final %>%
  filter_all(any_vars(str_detect(., "Neurocognitive Disorders"))) 
bipolar <- cc_final %>%
  filter_all(any_vars(str_detect(., "Bipolar and Related Disorders"))) 
eating_dis <- cc_final %>%
  filter_all(any_vars(str_detect(., "Feeding and Eating Disorders")))
# collapse into a list
dx_list <- list(adhd, asd, anxiety, disruptive, depressive, learning, communication, trauma, 
                intellectual_dis, motor_dis, ocd, scz, neurocog, bipolar, eating_dis)
# turn into binary 0/1 for each dx
dx <- map_dfc(dx_list, function(x){
  ifelse(test = cc_final$Anonymized.ID %in% x$Anonymized.ID, 
         yes = 1, no = 0)
}) %>%
  set_colnames(c("adhd", "asd", "anxiety", "disruptive", "depressive", "learning", "communication", 
                 "trauma", "intellectual_dis", "motor_dis", "ocd", "scz", "neurocog", 
                 "bipolar", "eating_dis")) %>%
  mutate(Anonymized.ID = cc_final$Anonymized.ID) %>%
  select(Anonymized.ID, everything())

##### once i figure out which participants to include, run code below...
# dx is a df of participants x dx categories, with binary 0/1 for absence/presence of dx
# split dx by tertile
bottom_dx <- dx[dx$Anonymized.ID %in% bottom33$Anonymized.ID, ]
top_dx <- dx[dx$Anonymized.ID %in% top33$Anonymized.ID, ]

# put into 15 2x2 matrices of high/low IQ x presence/absence of Dx **** corrected p = 0.00333
# ADHD !!
chisq.test(as.matrix(cbind(highIQ = table(top_dx$adhd), lowIQ = table(bottom_dx$adhd)))) #NS, surprisingly
# ASD !!
chisq.test(as.matrix(cbind(highIQ = table(top_dx$asd), lowIQ = table(bottom_dx$asd)))) # sig more in low IQ
# anxiety !!
chisq.test(as.matrix(cbind(highIQ = table(top_dx$anxiety), lowIQ = table(bottom_dx$anxiety)))) # sig more in high IQ
# depressive
chisq.test(as.matrix(cbind(highIQ = table(top_dx$depressive), lowIQ = table(bottom_dx$depressive)))) # sig more in high IQ
# disruptive !!
chisq.test(as.matrix(cbind(highIQ = table(top_dx$disruptive), lowIQ = table(bottom_dx$disruptive)))) # sig more in high
# learning !!
chisq.test(as.matrix(cbind(highIQ = table(top_dx$learning), lowIQ = table(bottom_dx$learning)))) # sig more in low
# communication !!
chisq.test(as.matrix(cbind(highIQ = table(top_dx$communication), lowIQ = table(bottom_dx$communication)))) # sig more in low
#intellectual disability !!
# bc there were no partic with ID in high IQ, need to manually enter otherwise chisq is wrong
chisq.test(as.matrix(cbind(highIQ = c(865,0), lowIQ = table(bottom_dx$intellectual_dis)))) # sig more in low
# OCD
chisq.test(as.matrix(cbind(highIQ = table(top_dx$ocd), lowIQ = table(bottom_dx$ocd)))) #NS
# trauma
chisq.test(as.matrix(cbind(highIQ = table(top_dx$trauma), lowIQ = table(bottom_dx$trauma)))) # NS
# motor
chisq.test(as.matrix(cbind(highIQ = table(top_dx$motor_dis), lowIQ = table(bottom_dx$motor_dis)))) # NS

### NOTE: the following Dx categories had too few participants with that Dx to assess:
# bipolar
# SCZ
# eating disorder
# neurocog




################################################
### SES MODERATION OF RACE X IQ RELATIONSHIP ### 
################################################
# regress SES out of IQ
linmodel <- lm(WISC_FSIQ ~ as.numeric(Barratt_Total), data = t1data)
# now, make a df of all the subscale scores that only include the observations that were included in the lm (some were left out bc of incomplete data)
data_ses_reg <- t1data[rownames(linmodel[["model"]]), ] %>%
  # add the residualized IQ score column
  mutate(fsiq_resid = linmodel[["residuals"]])
# check entered data in correct order; if not true stop the code
assertthat::assert_that(all(linmodel[["model"]]["WISC_FSIQ"] == data_ses_reg$WISC_FSIQ))
# remove participants who self-report American Indian/Alskan Native or Native Hawaiian, 
# bc 0 cases in the top IQ group, or who don't have race data
data_ses_reg <- data_ses_reg %>%
  filter(Child_Race != "American Indian/Alaskan Native" & Child_Race != "Native Hawaiian/Other Pacific Islander", 
         Child_Race != "Indian" & Child_Race != "Native American Indian" & Child_Race != "Unknown"
         & Child_Race != "Not Reported")
# split into tertiles based on SES-regressed IQ data
top33_ses <- data_ses_reg %>%
  slice_max(., fsiq_resid, prop = 0.33)
bottom33_ses <- data_ses_reg %>%
  slice_min(., fsiq_resid, prop = 0.33)
# chi-sq test of distribution of races in high and low iq groups after controlling for SES
chisq.test(as.matrix(cbind(highIQ = table(top33_ses$Child_Race), 
                           lowIQ = table(bottom33_ses$Child_Race)))) # p < 0.001

# compare to without controlling for SES
t1data_clean <- t1data %>%
  filter(Child_Race != "American Indian/Alaskan Native" & Child_Race != "Native Hawaiian/Other Pacific Islander", 
         Child_Race != "Indian" & Child_Race != "Native American Indian" & Child_Race != "Unknown"
         & Child_Race != "Not Reported")
( chisqres <- chisq.test(as.matrix(cbind(highIQ = table(t1data_clean$Child_Race[t1data_clean$group == "highIQ"]), 
                                         lowIQ = table(t1data_clean$Child_Race[t1data_clean$group == "lowIQ"])))) )
# residuals 
chisqres[["stdres"]]
# count matrix
ctmat <- as.matrix(cbind(highIQ = table(t1data_clean$Child_Race[t1data_clean$group == "highIQ"]), 
                lowIQ = table(t1data_clean$Child_Race[t1data_clean$group == "lowIQ"])))
# find percentages
percmat <- data.frame(highIQ = (ctmat[,1]/nrow(t1data_clean[t1data_clean$group == "highIQ",])*100), 
                      lowIQ = (ctmat[,2]/nrow(t1data_clean[t1data_clean$group == "lowIQ",])*100))

############################
### WSBM OF 55 IQ GROUPS ### 
############################
library(blockmodels)

# split into groups of 50 partic by iq
network_50 <- subscales %>%
  mutate(grp = as.factor(ntile(WISC_FSIQ, 55))) %>%
  group_split(grp)

# run wsbm on each group of 50, setting optimal number of clusters =9
WSBM_FXN_9 <- function(data){ # input x is df of participants x subscales
  # calculate correlation matrix 
  corrmat <- data %>%
    select(CBCL_AD_T:ICU_SR_Tot) %>%
    cor(., use = "pairwise.complete.obs") 
  # replace Inf values with 999, otherwise WSBM can't run
  corrmat[corrmat == Inf] <- 999
  # run WSBM 
  sbm <- BM_gaussian(membership_type = "SBM_sym", adj = corrmat)
  sbm$estimate()
  # determine optimum number of clusters // changed this to be all = 9, the modal number
  opt <- 9
  # get cluster assignments of each scale
  cluster_assgn <- data.frame(apply(sbm$memberships[[opt]]$Z, 1, function(x) which.max(x)))
  # calc avg iq score and avg corr and then put into df with subscale assignments
  output <- data.frame(avg_iq = mean(data$WISC_FSIQ, na.rm = TRUE),
                   avg_corr = mean(corrmat[upper.tri(corrmat, diag = FALSE)], na.rm = TRUE), 
                   number_blocks = max(cluster_assgn),
                   t(cluster_assgn)) %>%
    set_rownames("grp") %>%
    set_colnames(c("avg_iq", "avg_corr", "number_blocks", colnames(corrmat)))
  return(output)
}

# map fxn over each iq group to get df of iq group x vars (avg iq, avg corr, subscale assignments)
wsbm_by_iq <- map_dfr(network_50[1:55], function(d) WSBM_FXN_9(d))

# plot mean iq x correlation 
cor.test(wsbm_by_iq$avg_corr, wsbm_by_iq$avg_iq)
ggplot(wsbm_by_iq, aes(avg_iq, avg_corr)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Group average IQ", y = "Average pairwise correlation") +
  scale_x_continuous(n.breaks = 10) +
  scale_y_continuous(n.breaks = 6)
  # plot mean IQ x number of wsbm blocks
ggplot(wsbm_by_iq, aes(avg_iq, number_blocks)) +
  geom_jitter(width = 0.1, height = 0.1) + # add some jitter to be able to see each distinct point
  labs(x = "Group average IQ", y = "Optimal number of WSBM clusters") +
  scale_x_continuous(n.breaks = 10) 

# Hungarian method for determining maximum overlap of blocks

##### THIS CODE TAKEN FROM https://www.r-bloggers.com/2012/11/matching-clustering-solutions-using-the-hungarian-method/ 
##### ALL CREDIT GOES TO ORIGINAL AUTHOR, ROBERTO ROSLER

require(spdep)
require(rgdal)
require(maptools)
require(ggplot2)
require(plyr)
library(grid)

# function for determining best matches between clustering solutions in each group of 50 partic
# labels from cluster A will be matched on the labels from cluster B
minWeightBipartiteMatching <- function(clusteringA, clusteringB) {
  require(clue)
  idsA <- unique(clusteringA)  # distinct cluster ids in a
  idsB <- unique(clusteringB)  # distinct cluster ids in b
  nA <- length(clusteringA)  # number of instances in a
  nB <- length(clusteringB)  # number of instances in b
  if (length(idsA) != length(idsB) || nA != nB) {
    stop("number of cluster or number of instances do not match")
  }
  nC <- length(idsA)
  tupel <- c(1:nA)
  # computing the distance matrix
  assignmentMatrix <- matrix(rep(-1, nC * nC), nrow = nC)
  for (i in 1:nC) {
    tupelClusterI <- tupel[clusteringA == i]
    solRowI <- sapply(1:nC, function(i, clusterIDsB, tupelA_I) {
      nA_I <- length(tupelA_I)  # number of elements in cluster I
      tupelB_I <- tupel[clusterIDsB == i]
      nB_I <- length(tupelB_I)
      nTupelIntersect <- length(intersect(tupelA_I, tupelB_I))
      return((nA_I - nTupelIntersect) + (nB_I - nTupelIntersect))
    }, clusteringB, tupelClusterI)
    assignmentMatrix[i, ] <- solRowI
  }
  # optimization
  result <- solve_LSAP(assignmentMatrix, maximum = FALSE)
  attr(result, "assignmentMatrix") <- assignmentMatrix
  return(result)
}


# function that will perform the hungarian algorithm on each group of wsbm solutions
## NOTE that input d, all clustering solutions must contain the same number of clusters
wsbm_by_iq_list <- wsbm_by_iq %>%
  group_split(number_blocks) 
wsbm_by_iq_list <- map(wsbm_by_iq_list, function(d) d %>% select(!c(avg_iq, avg_corr, number_blocks)))
# function
HUNGARIAN_FXN <- function(d, matchon){ # input d is a df of iq groups x subscales in which each cell = cluster assignment of that subscale; matchon is seed row for matching
  output <- data.frame(matrix(NA, ncol = 73)) %>%
    set_colnames(colnames(d))
  for (i in 1:nrow(d)){
    clusterA <- as.numeric(d[i, 4:73]) # iterate across each row...
    names(clusterA) <- colnames(d)[4:73]
    clusterA_orig <- clusterA # need another structure with the results to use in the tmp loop to avoid overwriting
    clusterB <- as.numeric(d[matchon, 4:73]) # match clusters to solution in seed row
    names(clusterB) <- colnames(d)[4:73]
    # perform hungarian algorithm
    matching <- minWeightBipartiteMatching(clusterA, clusterB) 
    # change clusterA numbers to be the same as clusterB
    tmp <- sapply(1:length(matching), function(x) {
      clusterA[which(clusterA_orig == x)] <<- matching[x]
    })
    # output clusterA, and will eventually write into a df exactly like wsbm_by_iq, but now with matching # assignments
    output[i, 4:73] <- clusterA
  }
  # add avg iq, avg corr, and number of blocks to output df
  output[, 1:3] <- d[, 1:3]
  # return the output
  return(output)
}

# applying to 9-cluster solutions, matching on middle row
wsbm_by_iq_final <- HUNGARIAN_FXN(wsbm_by_iq, matchon = 27)
# plot
colors = structure(brewer.pal(n = 10, name = "Spectral"),
                   names = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")) # put 0 at the front so that the darkest red isn't used
Heatmap(as.matrix(wsbm_by_iq_final[, 4:73]), cluster_rows = FALSE, col = colors, 
        name = "Cluster", row_names_gp = gpar(fontsize = 9))






