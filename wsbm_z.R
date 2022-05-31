#####################################
#####################################
##### WSBM OF DIFFERENCE MATRIX #####
#####################################
#####################################

### SET UP 
library(dplyr)
library(tidyverse)
library(blockmodels)
library(superheat)
library(ComplexHeatmap)
setwd("~/Desktop/HBN")
theme_set(theme_bw())

### LOAD IN DIFFERENCE MATRIX --ALL WITH Z-TRANSFORMED NUMBERS
diffmat <- read.csv("./modified_data/diff_mat_WISC_FSIQ_33_ztransform.csv") %>%
  column_to_rownames("X") %>%
  as.matrix(.)
observed <- read.csv("./modified_data/deltar_percentiles_WISC_FSIQ_33_ztransform.csv")
# and null distribution
nulldist <- read.csv("./modified_data/null_delta_r_distribution_33_ztransform.csv")
# and subscales and human-readable names
subscales <- read.csv("./modified_data/ninth_release_iq_analyses2.csv")
new_names_heatmap <- read.csv("./modified_data/humanreadablenames.csv")
new_names_heatmap <- structure(new_names_heatmap$x, names = new_names_heatmap$X)
# add color and instrument assignments
inst_colors_df <- read.csv("./modified_data/instrument_color_assignments_df.csv")
inst_colors_vec <- read.csv("./modified_data/instrument_color_assignments70.csv")
inst_colors_vec <- structure(inst_colors_vec$x, names = inst_colors_vec$X)

#### RUN WSBM
diffmat_sbm <- BM_gaussian(membership_type = "SBM_sym", adj = diffmat)
diffmat_sbm$estimate()
# optimum number of clusters 
which.max(diffmat_sbm$ICL) # answer = 9
# cluster assignments 
diffmat_clust <- apply(diffmat_sbm$memberships[[9]]$Z, 1, function(x) which.max(x))
# determine block order by mean delta-r value of block 
diffmat_block_reorder <- c(1:9)[order(colMeans(diffmat_sbm$model_parameters[[9]]$mu))]
# new membership vector
diffmat_clust_reordered <- as.factor(diffmat_clust)
levels(diffmat_clust_reordered) <- rank(colMeans(diffmat_sbm$model_parameters[[9]]$m))
diffmat_clust_reordered <- as.numeric(as.character(diffmat_clust_reordered))
# df with cluster assignment of each variable
diffmat_clus_labs <- data.frame(scale = colnames(diffmat), sbm_cluster = diffmat_clust, 
                                sbm_cluster_ordered = diffmat_clust_reordered) %>%
  arrange(diffmat_clust_reordered)
# names for each block
diffmat_clus_labs$names <- factor(diffmat_clus_labs$sbm_cluster_ordered)
levels(diffmat_clus_labs$names) <- c("Aggressive behavior and RRBIs", 
                                     "SR Internalizing", "SR Aggressive and externalizing",
                                     "Total psychopathology, social, externalizing", 
                                     "SR Hyperactive/inattentive", "PR Internalizing", 
                                     "SR Anxiety", "Psychopathy", "Repetitive behaviors")
arrange(diffmat_clus_labs, sbm_cluster_ordered)
# visualize SBM
superheat(diffmat, membership.rows = diffmat_clust_reordered, membership.cols = diffmat_clust_reordered, 
          title = 'WSBM clustering of tertile difference matrix', 
          bottom.label = "variable", bottom.label.text.size = 2, 
          left.label.text.size = 2, bottom.label.text.angle = 90, 
          heat.pal = c("#276419", "white", "#8E0152")) # assign colors to values

### SAME THING USING COMPLEX HEATMAPS SYNTAX 
fixedcolors <- circlize::colorRamp2(breaks = c(0.3, 0, -0.3), colors = c("#276419", "white", "#8E0152"))
col2 <- circlize::colorRamp2(breaks = c(1, 0, -1), colors = c("#276419", "white", "#8E0152"))
nn <- read.csv("./modified_data/humanreadablenames.csv")
new_names_heatmap <- structure(nn$x, names = nn$X)
Heatmap(diffmat, name = "WSBM", col = fixedcolors, 
        # assign to WSBM block clusters
        row_split = diffmat_clust_reordered, column_split = diffmat_clust_reordered,
        row_title = NULL, column_title = NULL,
        # don't show dendrograms
        show_row_dend = FALSE, show_column_dend = FALSE, show_heatmap_legend = FALSE,
        # don't cluster slices, so that they show up in order of mean correlation with other blocks as in WSBM
        cluster_row_slices = FALSE, cluster_column_slices = FALSE, cluster_rows = FALSE, cluster_columns = FALSE,
        # human-readable names
        row_labels = new_names_heatmap[colnames(diffmat)], column_labels = new_names_heatmap[colnames(diffmat)],
        # square
        height = unit(16.5, "cm"), width = unit(16.5, "cm"),
        row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 8)) 

# add instrument and human-readable names and save WSBM groupings
diffmat_clus_labs <- diffmat_clus_labs %>%
  rowwise() %>%
  mutate(instrument = inst_colors_df$instrument[inst_colors_df$subscale %in% scale], 
         hr_name = inst_colors_df$hr_name[inst_colors_df$subscale %in% scale])
write.csv(diffmat_clus_labs, file = "./modified_data/WSBM_groupings.csv") # Table S1 in paper
           
           
### NULL DISTRIBUTION FOR WSBM BLOCKS

# df tracking which edges are in which block
edge_blocks <- data.frame(edge = nulldist$X) %>%
  # break edge into each subscale
  separate(., col = edge, into = c("subscale1", "subscale2"), sep = "_x_", remove = FALSE) %>%
  # add var for which block that subscale is in
  rowwise() %>%
  mutate(subscale1_block = diffmat_clus_labs$sbm_cluster_ordered[which(diffmat_clus_labs$scale %in% subscale1)], 
         subscale2_block = diffmat_clus_labs$sbm_cluster_ordered[which(diffmat_clus_labs$scale %in% subscale2)]) %>%
  ungroup() %>%
  select(edge, starts_with("subscale1"), starts_with("subscale2")) %>%
  # combine each edge's 2 blocks to get "edge blocks" 
  mutate(edge_block = ifelse(test = subscale1_block <= subscale2_block, 
                             yes = paste0(subscale1_block, subscale2_block), 
                             no = paste0(subscale2_block, subscale1_block))) %>%
  arrange(edge_block)
write.csv(edge_blocks, file = "./modified_data/WSBM_groupings_edges.csv")
# add edge block column to observed and null dist dfs
observed <- observed %>% 
  rowwise() %>%
  mutate(edge_block = edge_blocks$edge_block[which(edge_blocks$edge %in% pair)]) %>%
  select(X, edge_block, everything()) %>%
  ungroup()
nulldist <- nulldist %>%
  rowwise() %>%
  mutate(edge_block = edge_blocks$edge_block[which(edge_blocks$edge %in% X)]) %>%
  select(X, edge_block, everything()) %>%
  ungroup()

# df with block average delta r for each iteration of the null dist; rows are edge blocks and cols are each of the 10k iterations
block_nulls <- nulldist %>%
  group_by(edge_block) %>%
  summarise_if(.tbl = ., .predicate = is.numeric, mean) %>%
  ungroup()
# and for observed values
block_observed <- observed %>%
  group_by(edge_block) %>%
  summarise(block_avg = mean(observed_diff), 
            number_edges = n()) %>%
  ungroup()

# calculate percentile at which observed block avg falls compared to wsbm null dist
block_nulls_copy <- block_nulls %>%
  column_to_rownames("edge_block")
assertthat::assert_that(identical(rownames(block_nulls_copy), block_observed$edge_block))
percentiles_block <- data.frame(edge_block = matrix(rep(NA, nrow(block_nulls))), percentile = matrix(rep(NA, nrow(block_nulls))))
for (i in 1:nrow(block_nulls_copy)){
  percentiles_block$percentile[i] <- ecdf(block_nulls_copy[i, ])(block_observed$block_avg[i])
  percentiles_block$edge_block[i] <- block_observed$edge_block[i]
}

# add column for more descriptive block name
block_nulls <- block_nulls %>%
  # separate into 2 groups rather than 1 block
  mutate(edge_block = as.numeric(edge_block), 
         group1 = edge_block %/% 10, 
         group2 = edge_block %% 10) %>%
  # give descriptive group name, rather than number
  rowwise() %>%
  mutate(group1_name = unique(diffmat_clus_labs$names)[which(unique(diffmat_clus_labs$sbm_cluster_ordered) %in% group1)], 
         group2_name = unique(diffmat_clus_labs$names)[which(unique(diffmat_clus_labs$sbm_cluster_ordered) %in% group2)]) %>%
  ungroup() %>%
  # combine groups into blocks (2D)
  unite(col = "block_name", c(group1_name, group2_name), sep = "_x_", remove = TRUE) %>%
  select(!group1 & !group2) %>%
  select(edge_block, block_name, everything())

# add names to percentiles and observed dfs, too
percentiles_block <- percentiles_block %>%
  rowwise() %>%
  mutate(block_name = block_nulls$block_name[which(block_nulls$edge_block %in% edge_block)]) %>%
  ungroup() %>%
  # also add p-val and FDR-corrected p-value
  mutate(p_value = ifelse(test = percentile > 0.5, yes = 1 - percentile, no = percentile),
         adjusted_p_value = p.adjust(p_value, method = "fdr"))

block_observed <- block_observed %>%
  rowwise() %>%
  mutate(block_name = block_nulls$block_name[which(block_nulls$edge_block %in% edge_block)]) %>%
  ungroup()

# for Table S3
ts3 <- percentiles_block %>%
  rowwise() %>%
  mutate(block_avg_deltaz = block_observed$block_avg[block_observed$edge_block %in% edge_block]) %>%
  ungroup() %>%
  select(block_name, block_avg_deltaz, p_value) %>%
  mutate(block_name = str_replace(string = block_name, pattern = "_x_", replacement = " x "), 
         block_avg_deltaz = round(block_avg_deltaz, 3), 
         p_value = round(p_value, 3))
# save
write.csv(percentiles_block, file = "./modified_data/WSBM_percentiles_ztransform.csv")
write.csv(block_observed, file = "./modified_data/WSBM_block_avgs_ztransform.csv")
write.csv(ts3, file = "./figures/tablesupp3.csv")
which(percentiles_block$p_value <= 0.025) %>% length() # 14 blocks sig alpha =0.05
which(percentiles_block$p_value <= (0.025/45)) %>% length()     # 3 blocks survive Bonferroni 

### VISUALIZE as facetted histograms 
full_join(block_nulls, percentiles_block[, 2:5], by = "block_name") %>%
  # add in observed values 
  rowwise() %>%
  mutate(obs = block_observed$block_avg[which(block_observed$edge_block %in% edge_block)]) %>%
  ungroup() %>%
  select(edge_block, obs, percentile, p_value, adjusted_p_value, everything()) %>%
  # long form for plotting: want each block and its observed avg delta r repeated 1,000 times, one obs for each iteration
  pivot_longer(cols = !c(edge_block, block_name, obs, percentile, p_value, adjusted_p_value), values_to = "avg_deltar") %>%
  # remove name, just gives the index of the randomized null corr
  select(!name) %>%
  # add a variable that shows whether the null value is greater than the observed value or not, to change alpha in plot
  mutate(sig = as.factor(ifelse(adjusted_p_value <= 0.025, TRUE, FALSE)), 
         greater = case_when(
           obs < 0 ~ ifelse(test = (avg_deltar <= obs), yes = TRUE, no = FALSE), # if obs is negative, shade down from obs value
           obs > 0 ~ ifelse(test = (avg_deltar >= obs), yes = TRUE, no = FALSE)  # if obs is positive, shade up from obs value
         )) %>%
  # plot -- changing alpha to show where obs value falls
  ggplot(., aes(avg_deltar)) +
  geom_histogram(bins = 100, 
                 # change color based on whether or not obs is sig; change transparency based on whether null is > or < than obs
                 aes(fill = sig, alpha = greater), show.legend = FALSE) +      
  # separate plot for each pair of subscales
  facet_wrap(~ block_name) +
  theme(strip.text.x = element_text(size = 6)) +
  # change colors to make significant ones stand out more
  scale_fill_manual(values = c("grey", "red")) +
  # change the strength of the transparency to be more visible, and then change legend labels
  scale_alpha_discrete(range = c(0.35, 0.9)) +
  # add in percentile and observed correlation to each facet as text
  geom_text(data = percentiles_block, aes(x = -0.11, y = 900, label = paste0("%ile = ", percentile)), size = 3) +
  geom_text(data = block_observed, aes(x = -0.10, y = 700, label = paste0("obs. = ", round(block_avg, 4))), size = 3)


### WSBM DIFF MAT SHOWING GREY FOR NS EDGE BLOCKS
# matrix of proper row assignments into WSBM blocks
row_groups <- matrix(data = as.numeric(diffmat_clust_reordered) *10, 
                     nrow = 70, ncol = 70, byrow = FALSE, 
                     dimnames = list(rownames(diffmat), colnames(diffmat)))
# matrix of proper column assignments into WSBM blocks
col_groups <- matrix(data = as.numeric(diffmat_clust_reordered), 
                     nrow = 70, ncol = 70, byrow = TRUE, 
                     dimnames = list(rownames(diffmat), colnames(diffmat)))
# adding matrices together will give correct edge block assignments to each edge
diffmat_block_assignments <- matrix(data = (row_groups + col_groups), 
                                    nrow = 70, ncol = 70, byrow = TRUE, 
                                    dimnames = list(rownames(diffmat), colnames(diffmat)))
# replace the delta r values with NA for all blocks in which the block avg delta r is NS
# vector with edge block numbers of NS blocks
ns_blocks <- as.vector(as.numeric(percentiles_block$edge_block[percentiles_block$p_value > 0.025]))
diffmat_ns <- diffmat
diffmat_ns[diffmat_block_assignments %in% ns_blocks] <- NA

### VISUALIZE WITH NS PAIRS GREYED OUT
fixedcolors <- circlize::colorRamp2(breaks = c(0.3, 0, -0.3), colors = c("#276419", "white", "#8E0152"))
( hm_wsbm <- Heatmap(diffmat_ns, name = "WSBM", 
        col = fixedcolors, 
        # assign to WSBM block clusters
        row_split = diffmat_clust_reordered, column_split = diffmat_clust_reordered,
        row_title = NULL, column_title = NULL,
        # don't show dendrograms
        show_row_dend = FALSE, show_column_dend = FALSE, show_heatmap_legend = FALSE,
        # don't cluster slices, so that they show up in order of mean correlation with other blocks as in WSBM
        cluster_row_slices = FALSE, cluster_column_slices = FALSE, cluster_rows = FALSE, cluster_columns = FALSE,
        # human-readable names
        row_labels = new_names_heatmap[colnames(diffmat_ns)], column_labels = new_names_heatmap[colnames(diffmat_ns)],
        # square
        height = unit(12.5, "cm"), width = unit(12.5, "cm"),
        row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 8)) )
# now with instrument colors
inst_legend_wsbm <- Heatmap(as.matrix(colnames(diffmat_ns)), name = "Instrument Name", 
                        col = inst_colors_vec, cluster_rows = FALSE, 
                        show_heatmap_legend = FALSE, height = unit(12.5, "cm"), width = unit(0.5, "cm"))
# add color blocks for WSBM groupings
# instrument color assignments
block_colors <- structure(brewer.pal(n = 9, name = "Set1"), 
                          names = c("1", "2", "3", "4", "5", "6", "7", "8", "9")) 
block_colors_df <- diffmat_clus_labs %>%
  select(scale, sbm_cluster_ordered, names, hr_name) %>%
  mutate(color = block_colors[sbm_cluster_ordered])
block_colors_vec <- structure(block_colors_df$color, names = block_colors_df$scale)
group_legend <- Heatmap(as.matrix(colnames(diffmat_ns)), name = "Block Name", 
                            col = block_colors_vec, cluster_rows = FALSE, 
                            show_heatmap_legend = FALSE, height = unit(12.5, "cm"), width = unit(0.5, "cm"))
figwsbm <- hm_wsbm + group_legend + inst_legend_wsbm + 
  # add in names bc adding 2nd heatmap deletes them
  rowAnnotation(rn = anno_text(new_names_heatmap[colnames(diffmat_ns)]))
draw(figwsbm)

# legend for Fig 3
block_colors_df %>%
  mutate(sbm_cluster_ordered = as.character(sbm_cluster_ordered)) %>%
ggplot(., aes(x = sbm_cluster_ordered, y = 0.1)) + 
  geom_point(aes(color = sbm_cluster_ordered), shape = 15, size = 3) + 
  theme_classic() + 
  theme(axis.title = element_blank(), axis.line = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), 
        plot.margin = unit(c(0,0,0,0), "cm")) +
  labs(color = "Block Name") +
  scale_color_manual(values = block_colors, labels = c("Aggressive behavior and RRBIs", 
                                                       "SR Internalizing", "SR Aggressive and externalizing",
                                                       "Total psychopathology, social, externalizing", 
                                                       "SR Hyperactive/inattentive", "PR Internalizing", 
                                                       "SR Anxiety", "Psychopathy", "Repetitive behaviors"))

### FACETTED SCATTERPLOTS TO SHOW CORR FOR BLOCK AVERAGES IN HIGH VS LOW IQ GROUPS 
# df of z-scored scores 
zscores <- apply(subscales[, 17:86], 2, function(x) (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)) %>%
  as.data.frame()
# add in IQ column and anon ID
zscores2 <- zscores %>%
  mutate(IQ = subscales$WISC_FSIQ, 
         ID = subscales$Anonymized.ID) %>%
  # add column to track what IQ tertile it's in
  mutate(tertile = case_when(
    IQ >= 107 ~ "High IQ", 
    IQ <= 89 ~ "Low IQ"
  )) %>%
  # put into long form, then group the subscales into WSBM blocks 
  pivot_longer(cols = !c(IQ, ID, tertile), names_to = "subscale", values_to = "zscore") %>%
  rowwise() %>%
  mutate(wsbm_block = diffmat_clus_labs$names[which(diffmat_clus_labs$scale %in% subscale)]) %>%
  ungroup() %>%
  # group by person and by wsbm block to calculate the average score for that block per person
  group_by(ID, wsbm_block) %>%
  mutate(avg = mean(zscore, na.rm = TRUE)) %>%
  ungroup() %>%
  # select one row per combination of ID and wsbm_block
  distinct(., ID, wsbm_block, .keep_all = TRUE) %>%
  select(ID, IQ, tertile, wsbm_block, avg) %>%
  # pivot to wide form
  pivot_wider(names_from = wsbm_block, values_from = avg) %>%
  # don't care about obs not in one of the upper or lower tertiles, so can remove
  filter(!is.na(tertile))

# plots
ggplot(zscores2, aes(`Externalizing/ASD`, `Compulsions/self-injury/panic`)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "lm", aes(color = tertile)) +
  scale_color_manual(values=c("red", "blue")) +
  labs(x = "Externalizing/social impairments") +
  facet_wrap(~ tertile) +
  theme(axis.title = element_text(size = 15), strip.text = element_text(size = 13)) +
  stat_fit_glance(method = "lm",
                  method.args = list(formula = y ~ x), 
                  label.y = "top", 
                  label.x = "right",
                  mapping = aes(label = sprintf('r~"="~%.3f',
                                                sqrt(stat(r.squared)))), 
                  parse = TRUE, size = 7)


### SLOPEGRAPH / DUMBBELL PLOT ###
# pull out corrs in upper and lower tertiles
kept <- which(upper.tri(split.corr[["high"]], diag = FALSE), arr.ind = TRUE)   # indices of row and column names that were kept with upper.tri
tert_corrs <- data.frame(pair =  paste(dimnames(split.corr[["high"]])[[2]][kept[ , 2]], dimnames(split.corr[["high"]])[[1]][kept[ , 1]], sep = "_x_"), 
                         high_iq = split.corr[["high"]][upper.tri(split.corr[["high"]], diag = FALSE)], 
                         low_iq = split.corr[["low"]][upper.tri(split.corr[["low"]], diag = FALSE)]) %>%
  # only keep those edges in edge block 24
  filter(pair %in% edge_blocks$edge[edge_blocks$edge_block == 24]) %>%
  # long form
  pivot_longer(., cols = !pair, names_to = "iq", values_to = "corr")
# plot slopegraph
ggplot(tert_corrs, aes(x = iq, y = corr, group = pair)) +
  geom_point() + 
  geom_line() +
  labs(x = "IQ Tertile", y = "Correlation") +
  scale_x_discrete(labels = c("Upper", "Lower"))
# plot dumbbell plot
ggplot(tert_corrs, aes(x = corr, y = fct_reorder(pair, corr))) +
  geom_point(aes(col = iq)) +
  geom_line() +
  labs(x = "Correlation", y = "Edge") +
  scale_color_discrete(name = "IQ Tertile", labels = c("High IQ", "Low IQ")) +
  theme(axis.text.y = element_text(size = 10))



###################################
### WSBM USING IQ RESIDUAL DATA ###
###################################
### LOAD IN DIFFERENCE MATRIX 
diffmat_resid_iq <- read.csv("./modified_data/diff_mat_iq_regressed_z.csv") %>%
  column_to_rownames("X") %>%
  as.matrix(.)
# replace NA with 0 on diagonal blocks
diffmat_resid_iq[which(is.na(diffmat_resid_iq) == TRUE)] <- 0

#### RUN WSBM
diffmat_sbm_riq <- BM_gaussian(membership_type = "SBM_sym", adj = diffmat_resid_iq)
diffmat_sbm_riq$estimate()
# optimum number of clusters 
which.max(diffmat_sbm_riq$ICL) # answer = 8
# cluster assignments 
diffmat_clust_riq <- apply(diffmat_sbm_riq$memberships[[8]]$Z, 1, function(x) which.max(x))
# determine block order by mean delta-r value of block 
diffmat_block_reorder_riq <- c(1:8)[order(colMeans(diffmat_sbm_riq$model_parameters[[8]]$mu))]
# new membership vector
diffmat_clust_reordered_riq <- as.factor(diffmat_clust_riq)
levels(diffmat_clust_reordered_riq) <- rank(colMeans(diffmat_sbm_riq$model_parameters[[8]]$m))
diffmat_clust_reordered_riq <- as.numeric(as.character(diffmat_clust_reordered_riq))
# df with cluster assignment of each variable
diffmat_clus_labs_riq <- data.frame(scale = colnames(diffmat_resid_iq), sbm_cluster = diffmat_clust_riq, 
                                sbm_cluster_ordered = diffmat_clust_reordered_riq) %>%
  arrange(diffmat_clust_reordered_riq)
# names for each block
diffmat_clus_labs_riq$names <- factor(diffmat_clus_labs_riq$sbm_cluster_ordered)
levels(diffmat_clus_labs_riq$names) <- c("Externalizing/aggression and social problems", "SR Int", "Misc1",
                                  "Anxiety/affective prob", "Hyperactive/inattentive and social prob", 
                                  "Misc2", "Psychopathy", "Repetitive behaviors")


Heatmap(diffmat_resid_iq, name = "WSBM", col = fixedcolors,
        # assign to WSBM block clusters
        row_split = diffmat_clust_reordered_riq, column_split = diffmat_clust_reordered_riq, 
        # don't show dendrograms
        show_row_dend = FALSE, show_column_dend = FALSE, 
        # don't cluster slices, so that they show up in order of mean correlation with other blocks as in WSBM
        cluster_row_slices = FALSE, cluster_column_slices = FALSE, 
        row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 7))



### NULL DISTRIBUTION FOR WSBM BLOCKS
# load in iq residualized null dist and percentiles
nulldist_riq <- read.csv("./modified_data/null_delta_r_distribution_iq_resid.csv")
observed_riq <- read.csv("./modified_data/deltar_percentiles_iq_resid.csv")
# df tracking which edges are in which block
edge_blocks_riq <- data.frame(edge = nulldist_riq$X) %>%
  # break edge into each subscale
  separate(., col = edge, into = c("subscale1", "subscale2"), sep = "_x_", remove = FALSE) %>%
  # add var for which block that subscale is in
  rowwise() %>%
  mutate(subscale1_block = diffmat_clus_labs_riq$sbm_cluster_ordered[which(diffmat_clus_labs_riq$scale %in% subscale1)], 
         subscale2_block = diffmat_clus_labs_riq$sbm_cluster_ordered[which(diffmat_clus_labs_riq$scale %in% subscale2)]) %>%
  ungroup() %>%
  select(edge, starts_with("subscale1"), starts_with("subscale2")) %>%
  # combine each edge's 2 blocks to get "edge blocks" 
  mutate(edge_block = ifelse(test = subscale1_block <= subscale2_block, 
                             yes = paste0(subscale1_block, subscale2_block), 
                             no = paste0(subscale2_block, subscale1_block))) %>%
  arrange(edge_block)

# add edge block column to observed and null dist dfs
observed_riq <- observed_riq %>% 
  rowwise() %>%
  mutate(edge_block = edge_blocks_riq$edge_block[which(edge_blocks_riq$edge %in% pair)]) %>%
  select(X, edge_block, everything()) %>%
  ungroup()
nulldist_riq <- nulldist_riq %>%
  rowwise() %>%
  mutate(edge_block = edge_blocks_riq$edge_block[which(edge_blocks_riq$edge %in% X)]) %>%
  select(X, edge_block, everything()) %>%
  ungroup()

# df with block average delta r for each iteration of the null dist; rows are edge blocks and cols are each of the 10k iterations
block_nulls_riq <- nulldist_riq %>%
  group_by(edge_block) %>%
  summarise_if(.tbl = ., .predicate = is.numeric, mean) %>%
  ungroup()
# and for observed values
block_observed_riq <- observed_riq %>%
  group_by(edge_block) %>%
  summarise(block_avg = mean(observed_diff), 
            number_edges = n()) %>%
  ungroup()

# calculate percentile at which observed block avg falls compared to wsbm null dist
block_nulls_copy_riq <- block_nulls_riq %>%
  column_to_rownames("edge_block")
assertthat::assert_that(identical(rownames(block_nulls_copy_riq), block_observed_riq$edge_block))
percentiles_block_riq <- data.frame(edge_block = matrix(rep(NA, nrow(block_nulls_riq))), percentile = matrix(rep(NA, nrow(block_nulls_riq))))
for (i in 1:nrow(block_nulls_copy_riq)){
  percentiles_block_riq$percentile[i] <- ecdf(block_nulls_copy_riq[i, ])(block_observed_riq$block_avg[i])
  percentiles_block_riq$edge_block[i] <- block_observed_riq$edge_block[i]
}


# add column for more descriptive block name
block_nulls_riq <- block_nulls_riq %>%
  # separate into 2 groups rather than 1 block
  mutate(edge_block = as.numeric(edge_block), 
         group1 = edge_block %/% 10, 
         group2 = edge_block %% 10) %>%
  # give descriptive group name, rather than number
  rowwise() %>%
  mutate(group1_name = unique(diffmat_clus_labs_riq$names)[which(unique(diffmat_clus_labs_riq$sbm_cluster_ordered) %in% group1)], 
         group2_name = unique(diffmat_clus_labs_riq$names)[which(unique(diffmat_clus_labs_riq$sbm_cluster_ordered) %in% group2)]) %>%
  ungroup() %>%
  # combine groups into blocks (2D)
  unite(col = "block_name", c(group1_name, group2_name), sep = "_x_", remove = TRUE) %>%
  select(!group1 & !group2) %>%
  select(edge_block, block_name, everything())

# add names to percentiles and observed dfs, too
percentiles_block_riq <- percentiles_block_riq %>%
  rowwise() %>%
  mutate(block_name = block_nulls_riq$block_name[which(block_nulls_riq$edge_block %in% edge_block)]) %>%
  ungroup() %>%
  # also add p-val 
  mutate(p_value = ifelse(test = percentile > 0.5, yes = 1 - percentile, no = percentile))

block_observed_riq <- block_observed_riq %>%
  rowwise() %>%
  mutate(block_name = block_nulls_riq$block_name[which(block_nulls_riq$edge_block %in% edge_block)]) %>%
  ungroup()

# save
write.csv(percentiles_block_riq, file = "./modified_data/WSBM_percentiles_iq_resid.csv")
write.csv(block_observed_riq, file = "./modified_data/WSBM_block_avgs_iq_resid.csv")

which(percentiles_block_riq$p_value <= (0.025/36)) %>% length()     # 2 blocks survive Bonferroni 

### WSBM DIFF MAT SHOWING GREY FOR NS EDGE BLOCKS
# matrix of proper row assignments into WSBM blocks
row_groups_riq <- matrix(data = as.numeric(diffmat_clust_reordered_riq) *10, 
                     nrow = 70, ncol = 70, byrow = FALSE, 
                     dimnames = list(rownames(diffmat_resid_iq), colnames(diffmat_resid_iq)))
# matrix of proper column assignments into WSBM blocks
col_groups_riq <- matrix(data = as.numeric(diffmat_clust_reordered_riq), 
                     nrow = 70, ncol = 70, byrow = TRUE, 
                     dimnames = list(rownames(diffmat_resid_iq), colnames(diffmat_resid_iq)))
# adding matrices together will give correct edge block assignments to each edge
diffmat_block_assignments_riq <- matrix(data = (row_groups_riq + col_groups_riq), 
                                    nrow = 70, ncol = 70, byrow = TRUE, 
                                    dimnames = list(rownames(diffmat_resid_iq), colnames(diffmat_resid_iq)))
# replace the delta r values with NA for all blocks in which the block avg delta r is NS
# vector with edge block numbers of NS blocks
ns_blocks_riq <- as.vector(as.numeric(percentiles_block_riq$edge_block[percentiles_block_riq$p_value > 0.025]))
diffmat_ns_riq <- diffmat_resid_iq
diffmat_ns_riq[diffmat_block_assignments_riq %in% ns_blocks_riq] <- NA

### VISUALIZE WITH NS PAIRS GREYED OUT
fixedcolors <- circlize::colorRamp2(breaks = c(0.3, 0, -0.3), colors = c("#276419", "white", "#8E0152"))
( hm_ns_riq <- Heatmap(diffmat_ns_riq, name = "WSBM", 
        col = fixedcolors, row_title = NULL, column_title = NULL,
        # assign to WSBM block clusters
        row_split = diffmat_clust_reordered_riq, column_split = diffmat_clust_reordered_riq, 
        # don't show dendrograms
        show_row_dend = FALSE, show_column_dend = FALSE, show_heatmap_legend = FALSE,
        # don't cluster slices, so that they show up in order of mean correlation with other blocks as in WSBM
        cluster_row_slices = FALSE, cluster_column_slices = FALSE, 
        cluster_rows = FALSE, cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 8), 
        # human-readable names
        row_labels = new_names_heatmap[colnames(diffmat_ns_riq)], column_labels = new_names_heatmap[colnames(diffmat_ns_riq)],
        # square
        height = unit(16.5, "cm"), width = unit(16.5, "cm")) )
# now with instrument colors
inst_legend_riq <- Heatmap(as.matrix(colnames(diffmat_ns_riq)), name = "Instrument Name", 
                            col = inst_colors_vec, cluster_rows = FALSE, 
                            show_heatmap_legend = FALSE, height = unit(16.5, "cm"), width = unit(0.5, "cm"))
figwsbm_riq <- hm_ns_riq + inst_legend_riq + 
  # add in names bc adding 2nd heatmap deletes them
  rowAnnotation(rn = anno_text(new_names_heatmap[colnames(diffmat_ns_riq)]))
draw(figwsbm_riq)










#########################################
### WSBM USING AGE/SEX REGRESSED DATA ###
#########################################
### LOAD IN DIFFERENCE MATRIX 
diffmat_resid_bas <- read.csv("./modified_data/diff_mat_age_sex_regressed_z.csv") %>%
  column_to_rownames("X") %>%
  as.matrix(.)
# replace NA with 0 on diagonal blocks
diffmat_resid_bas[which(is.na(diffmat_resid_bas) == TRUE)] <- 0

#### RUN WSBM
diffmat_sbm_bas <- BM_gaussian(membership_type = "SBM_sym", adj = diffmat_resid_iq)
diffmat_sbm_bas$estimate()
# optimum number of clusters 
which.max(diffmat_sbm_bas$ICL) # answer = 8
# cluster assignments 
diffmat_clust_bas <- apply(diffmat_sbm_bas$memberships[[8]]$Z, 1, function(x) which.max(x))
# determine block order by mean delta-r value of block 
diffmat_block_reorder_bas <- c(1:8)[order(colMeans(diffmat_sbm_bas$model_parameters[[8]]$mu))]
# new membership vector
diffmat_clust_reordered_bas <- as.factor(diffmat_clust_bas)
levels(diffmat_clust_reordered_bas) <- rank(colMeans(diffmat_sbm_bas$model_parameters[[8]]$m))
diffmat_clust_reordered_bas <- as.numeric(as.character(diffmat_clust_reordered_bas))
# df with cluster assignment of each variable
diffmat_clus_labs_bas <- data.frame(scale = colnames(diffmat_resid_bas), sbm_cluster = diffmat_clust_bas, 
                                    sbm_cluster_ordered = diffmat_clust_reordered_bas) %>%
  arrange(diffmat_clust_reordered_bas)
# names for each block
diffmat_clus_labs_bas$names <- factor(diffmat_clus_labs_bas$sbm_cluster_ordered)
levels(diffmat_clus_labs_bas$names) <- c("Externalizing/aggressive", "SR Int", "Misc1",
                                   "Anxiety/affective", "SR Inattentive/hyperactive and social prob", 
                                   "Misc2", "Psychopathy", "Repetitive behaviors")

Heatmap(diffmat_resid_bas, name = "WSBM", col = fixedcolors,
        # assign to WSBM block clusters
        row_split = diffmat_clust_reordered_bas, column_split = diffmat_clust_reordered_bas, 
        # don't show dendrograms
        show_row_dend = FALSE, show_column_dend = FALSE, 
        # don't cluster slices, so that they show up in order of mean correlation with other blocks as in WSBM
        cluster_row_slices = FALSE, cluster_column_slices = FALSE, 
        row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 7))



### NULL DISTRIBUTION FOR WSBM BLOCKS
# load in iq residualized null dist and percentiles
nulldist_bas <- read.csv("./modified_data/null_delta_r_distribution_agesex_resid.csv")
observed_bas <- read.csv("./modified_data/deltar_percentiles_agesex_resid.csv")
# df tracking which edges are in which block
edge_blocks_bas <- data.frame(edge = nulldist_bas$X) %>%
  # break edge into each subscale
  separate(., col = edge, into = c("subscale1", "subscale2"), sep = "_x_", remove = FALSE) %>%
  # add var for which block that subscale is in
  rowwise() %>%
  mutate(subscale1_block = diffmat_clus_labs_bas$sbm_cluster_ordered[which(diffmat_clus_labs_bas$scale %in% subscale1)], 
         subscale2_block = diffmat_clus_labs_bas$sbm_cluster_ordered[which(diffmat_clus_labs_bas$scale %in% subscale2)]) %>%
  ungroup() %>%
  select(edge, starts_with("subscale1"), starts_with("subscale2")) %>%
  # combine each edge's 2 blocks to get "edge blocks" 
  mutate(edge_block = ifelse(test = subscale1_block <= subscale2_block, 
                             yes = paste0(subscale1_block, subscale2_block), 
                             no = paste0(subscale2_block, subscale1_block))) %>%
  arrange(edge_block)

# add edge block column to observed and null dist dfs
observed_bas <- observed_bas %>% 
  rowwise() %>%
  mutate(edge_block = edge_blocks_bas$edge_block[which(edge_blocks_bas$edge %in% pair)]) %>%
  select(X, edge_block, everything()) %>%
  ungroup()
nulldist_bas <- nulldist_bas %>%
  rowwise() %>%
  mutate(edge_block = edge_blocks_bas$edge_block[which(edge_blocks_bas$edge %in% X)]) %>%
  select(X, edge_block, everything()) %>%
  ungroup()

# df with block average delta r for each iteration of the null dist; rows are edge blocks and cols are each of the 10k iterations
block_nulls_bas <- nulldist_bas %>%
  group_by(edge_block) %>%
  summarise_if(.tbl = ., .predicate = is.numeric, mean) %>%
  ungroup()
# and for observed values
block_observed_bas <- observed_bas %>%
  group_by(edge_block) %>%
  summarise(block_avg = mean(observed_diff), 
            number_edges = n()) %>%
  ungroup()

# calculate percentile at which observed block avg falls compared to wsbm null dist
block_nulls_copy_bas <- block_nulls_bas %>%
  column_to_rownames("edge_block")
assertthat::assert_that(identical(rownames(block_nulls_copy_bas), block_observed_bas$edge_block))
percentiles_block_bas <- data.frame(edge_block = matrix(rep(NA, nrow(block_nulls_bas))), percentile = matrix(rep(NA, nrow(block_nulls_bas))))
for (i in 1:nrow(block_nulls_copy_bas)){
  percentiles_block_bas$percentile[i] <- ecdf(block_nulls_copy_bas[i, ])(block_observed_bas$block_avg[i])
  percentiles_block_bas$edge_block[i] <- block_observed_bas$edge_block[i]
}

# add column for more descriptive block name
block_nulls_bas <- block_nulls_bas %>%
  # separate into 2 groups rather than 1 block
  mutate(edge_block = as.numeric(edge_block), 
         group1 = edge_block %/% 10, 
         group2 = edge_block %% 10) %>%
  # give descriptive group name, rather than number
  rowwise() %>%
  mutate(group1_name = unique(diffmat_clus_labs_bas$names)[which(unique(diffmat_clus_labs_bas$sbm_cluster_ordered) %in% group1)], 
         group2_name = unique(diffmat_clus_labs_bas$names)[which(unique(diffmat_clus_labs_bas$sbm_cluster_ordered) %in% group2)]) %>%
  ungroup() %>%
  # combine groups into blocks (2D)
  unite(col = "block_name", c(group1_name, group2_name), sep = "_x_", remove = TRUE) %>%
  select(!group1 & !group2) %>%
  select(edge_block, block_name, everything())

# add names to percentiles and observed dfs, too
percentiles_block_bas <- percentiles_block_bas %>%
  rowwise() %>%
  mutate(block_name = block_nulls_bas$block_name[which(block_nulls_bas$edge_block %in% edge_block)]) %>%
  ungroup() %>%
  # also add p-val and FDR-corrected p-value
  mutate(p_value = ifelse(test = percentile > 0.5, yes = 1 - percentile, no = percentile))

block_observed_bas <- block_observed_bas %>%
  rowwise() %>%
  mutate(block_name = block_nulls_bas$block_name[which(block_nulls_bas$edge_block %in% edge_block)]) %>%
  ungroup()

# save
write.csv(percentiles_block_bas, file = "./modified_data/WSBM_percentiles_agesex_resid.csv")
write.csv(block_observed_bas, file = "./modified_data/WSBM_block_avgs_agesex_resid.csv")

which(percentiles_block_bas$p_value <= (0.025/36)) %>% length()     # 1 block survives Bonferroni 

### WSBM DIFF MAT SHOWING GREY FOR NS EDGE BLOCKS
# matrix of proper row assignments into WSBM blocks
row_groups_bas <- matrix(data = as.numeric(diffmat_clust_reordered_bas) *10, 
                         nrow = 70, ncol = 70, byrow = FALSE, 
                         dimnames = list(rownames(diffmat_resid_bas), colnames(diffmat_resid_bas)))
# matrix of proper column assignments into WSBM blocks
col_groups_bas <- matrix(data = as.numeric(diffmat_clust_reordered_bas), 
                         nrow = 70, ncol = 70, byrow = TRUE, 
                         dimnames = list(rownames(diffmat_resid_bas), colnames(diffmat_resid_bas)))
# adding matrices together will give correct edge block assignments to each edge
diffmat_block_assignments_bas <- matrix(data = (row_groups_bas + col_groups_bas), 
                                        nrow = 70, ncol = 70, byrow = TRUE, 
                                        dimnames = list(rownames(diffmat_resid_bas), colnames(diffmat_resid_bas)))
# replace the delta r values with NA for all blocks in which the block avg delta r is NS
# vector with edge block numbers of NS blocks
ns_blocks_bas <- as.vector(as.numeric(percentiles_block_bas$edge_block[percentiles_block_bas$p_value > 0.025]))
diffmat_ns_bas <- diffmat_resid_bas
diffmat_ns_bas[diffmat_block_assignments_bas %in% ns_blocks_bas] <- NA

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
new_names_heatmap_bas <- structure(new_names, names = colnames(diffmat_bas))
### VISUALIZE WITH NS PAIRS GREYED OUT
fixedcolors <- circlize::colorRamp2(breaks = c(0.3, 0, -0.3), colors = c("#276419", "white", "#8E0152"))
( hm_ns_bas <-Heatmap(diffmat_ns_bas, name = "WSBM", 
        col = fixedcolors, row_title = NULL, column_title = NULL,
        # assign to WSBM block clusters
        row_split = diffmat_clust_reordered_bas, column_split = diffmat_clust_reordered_bas, 
        # don't show dendrograms
        show_row_dend = FALSE, show_column_dend = FALSE, show_heatmap_legend = FALSE,
        # don't cluster slices, so that they show up in order of mean correlation with other blocks as in WSBM
        cluster_row_slices = FALSE, cluster_column_slices = FALSE, 
        cluster_rows = FALSE, cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 8), 
        # human-readable names
        row_labels = new_names_heatmap_bas[colnames(diffmat_ns_bas)], column_labels = new_names_heatmap_bas[colnames(diffmat_ns_bas)],
        # square
        height = unit(16.5, "cm"), width = unit(16.5, "cm")) )
# now with instrument colors
observed_means_bas <- observed_means_bas %>%
  mutate(instrument = c(rep("CBCL", 11), rep("YSR", 11), rep("CSR", 5), rep("SRS", 8), rep("SCARED_PR", 6), 
                        rep("SCARED_SR", 6), rep("RBS", 6), rep("SDQ", 9), rep("ICU_PR", 4), rep("ICU_SR",4)))
inst_colors_bas <- observed_means_bas %>%
  select(subscale, instrument) %>%
  mutate(color = instrument_colors[instrument], 
         hr_name = new_names_heatmap_bas)
inst_colors_vec_bas <- structure(inst_colors_bas$color, names = inst_colors_bas$subscale)
inst_legend_bas <- Heatmap(as.matrix(colnames(diffmat_ns_bas)), name = "Instrument Name", 
                           col = inst_colors_vec_bas, cluster_rows = FALSE, 
                           show_heatmap_legend = FALSE, height = unit(16.5, "cm"), width = unit(0.5, "cm"))
figwsbm_bas <- hm_ns_bas + inst_legend_bas + 
  # add in names bc adding 2nd heatmap deletes them
  rowAnnotation(rn = anno_text(new_names_heatmap_bas[colnames(diffmat_ns_bas)]))
draw(figwsbm_bas)

##########################################################
### CHI-SQUARE TEST FOR SUBSCALE ASSIGNMENTS TO BLOCKS ###
##########################################################
# put all in same order
d <- data.frame(blocks_main = diffmat_clus_labs[order(diffmat_clus_labs$scale), 3], # selecting 3rd col, sbm_cluster_ordered
                blocks_bas = diffmat_clus_labs_bas[order(diffmat_clus_labs_bas$scale), 3], 
                blocks_riq = diffmat_clus_labs_riq[order(diffmat_clus_labs_riq$scale), 3])

# comparing main analyses to when controlling for age and sex in behavioral subscales
as.matrix(cbind(table(d$blocks_main), table(d$blocks_bas)))
with(d, table(d$blocks_main, d$blocks_bas)) %>% chisq.test()
# comparing main analyses to when controlling for SES and race in IQ score
with(d, table(blocks_main, blocks_riq)) %>% chisq.test()







#################################################
#################################################
### ANALYSES IN RESPONSE TO REVIEWER COMMENTS ###
#################################################
#################################################

########################################################################
### DETERMINING OVERLAP BETWEEN WSBM SIG BLOCKS AND INDIVIDUAL EDGES ###
########################################################################
# make sure the following vars are loaded in: percentiles, percentiles_block, block_observed, diffmat_clus_labs, edge_blocks

# binarize whether each block and each edge is sig vs.non-sig
blocks <- percentiles_block %>%
  rowwise() %>%
  mutate(block_avg = block_observed$block_avg[block_observed$edge_block %in% edge_block]) %>%
  ungroup() %>%
  mutate(sig_pos = ifelse(test = p_value <= 0.025 & block_avg > 0, yes = 1, no = 0),
         sig_neg = ifelse(test = p_value <= 0.025 & block_avg < 0, yes = 1, no = 0))
edges <- percentiles %>%
  mutate(sig_pos = ifelse(test = p_value <= 0.025 & observed_diff > 0, yes = 1, no = 0),
         sig_neg = ifelse(test = p_value <= 0.025 & observed_diff < 0, yes = 1, no = 0)) %>%
  # add info on which block each edge is in
  rowwise() %>%
  mutate(block = edge_blocks$edge_block[edge_blocks$edge %in% pair]) %>%
  ungroup()
# calculate number of significant edges (+ or -) in each edge block
number_sig <- edges %>%
  group_by(block) %>%
  summarise(n_pos = sum(sig_pos),
            n_neg = sum(sig_neg)) %>%
  ungroup() %>%
  # and add column for whether the block is significant in WSBM analyses
  mutate(block_sig_pos = blocks$sig_pos,
         block_sig_neg = blocks$sig_neg, 
         # and how many edges in the block total
         n_edges = block_observed$number_edges) %>%
  # only care about those blocks that are significant
  filter(block_sig_pos | block_sig_neg == 1)

# draw n edges 1000 times (n = number of total edges in each significant block)
output <- matrix(nrow = 14, ncol = 1000) # make empty df that will contain results
for (n in 1:1000){ # first for-loop will loop through the 1000 permutations 
  for (i in 1:nrow(number_sig)){ # nested for-loop will loop through the 14 significant blocks in the number_sig df
    # take a random sample of n edges, where n = number of edges in block i
    random_sample <- slice_sample(edges, n = number_sig$n_edges[i])
    # count how many of these random edges are significantly neg or pos, depending on whether the block was neg or pos
    val <- ifelse(test = number_sig$block_sig_pos[i] == 1, # if the block is significantly positive...
                  yes = sum(random_sample$sig_pos),        # sum how many of the random edges are sig pos
                  no = sum(random_sample$sig_neg))         # otherwise, sum how many are sig neg
    output[i, n] <- val
  }
}

# determine at what percentile the observed value falls 
## first, need to collapse n_pos and n_sig into one column
number_sig$n_sig <- ifelse(test = number_sig$block_sig_pos == 1, # if the block is significantly positive...
                           yes = number_sig$n_pos,     # sum how many of the random edges are sig pos
                           no = number_sig$n_neg)
## make empty vector for results
observed_count_percentile <- vector(mode = "double", length = 14)
for (i in 1:14){
  observed_count_percentile[i] <- ecdf(output[i, 1:1000])(number_sig$n_sig[i])
}
observed_count_percentile # 12/14 significant blocks contain significantly more significant edges than expected by chance (at p < 0.025... 14/14 if p < 0.05)

#########################
### EXCLUDING IQ < 50 ### 
#########################

low <- bottom33 %>%
  filter(WISC_FSIQ >=50)
# collapse into a list
low_vs_upper <- list(high = top33, low = low)

# calculate mean score for each subscale 
means_50 <- map_dfc(low_vs_upper, function(d){
  d %>%
    # select subscales
    select(CBCL_AD_T : ICU_SR_Tot) %>%
    # compute mean
    apply(., 2, mean, na.rm = TRUE) %>%
    data.frame(.) %>%
    magrittr::set_rownames(colnames(d[17:86]))
}) %>%
  # save the means data frame in this order, so i can use it to compare against the null distribution later on
  set_colnames(c("high", "low_50")) %>%
  rownames_to_column("subscale") %>%
  # add a var for the difference between mean score in low vs high IQ groups
  mutate(difference = high - low_50) 


# calculate standard deviation of each subscale in the high and low IQ groups
stdev_50 <- map_dfc(low_vs_upper, function(d){
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
  set_colnames(c("high", "low_50")) %>%
  rownames_to_column("subscale") %>%
  # compute the pooled standard deviation, which will be used to caluclate cohen's d
  mutate(pooled_sd = sqrt(((high ^ 2) + (low_50 ^ 2)) / 2) )

# divide the difference in means by the pooled standard deviation
# first check that rows are in the same order
assertthat::assert_that(identical(means_50$subscale, stdev_50$subscale))
cohens_d_50 <- data.frame(subscale = means_50$subscale, d = means_50$difference / stdev_50$pooled_sd)

### COMPARE TO MAIN RESULTS ###
cohens_d <- read.csv("./modified_data/cohens_d_WISC_FSIQ_33.csv")
cor(cohens_d$d, cohens_d_50$d) # r=0.999 very high

# calculate correlation matrix for each tertile, then z-transform
split.corr_50 <- map(low_vs_upper, function(d){
  d %>% 
    select(CBCL_AD_T : ICU_SR_Tot) %>%
    cor(., use = "pairwise.complete.obs") %>% 
    # z-transform
    fisherz()
})

# subtract LOW from HIGH to make difference matrix
assertthat::assert_that(identical(colnames(split.corr_50[["low"]]),  colnames(split.corr_50[["high"]]))) # check row and col names are identical
assertthat::assert_that(identical(rownames(split.corr_50[["low"]]),  rownames(split.corr_50[["high"]])))
# make diagonal values 0 instead of Inf
split.corr_50[[1]][which(is.infinite(split.corr_50[[1]]) == TRUE)]<- 0
split.corr_50[[2]][which(is.infinite(split.corr_50[[2]]) == TRUE)]<- 0
# subtract
diffmat_50 <- split.corr_50[["high"]] - split.corr_50[["low"]] 

### COMPARE TO MAIN RESULTS ###
lmh_z <- read.csv("./modified_data/diff_mat_WISC_FSIQ_33_ztransform.csv") %>%
  column_to_rownames("X")
cor(diffmat_50[upper.tri(diffmat_50, dia = FALSE)], lmh_z[upper.tri(lmh_z, dia = FALSE)]) # r=0.9999 very high


