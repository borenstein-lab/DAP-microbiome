# Takes a feature table with a "sample_id" column and replaces sample_id's 
#  with dog_id's.
# In case the table is a distance matrix, i.e. includes sample_id's in the 
#  column names as well, these are replaced to.
map_sample_id_to_dap_dog_id <- function(feat_table, 
                                        sample_map, 
                                        rename_rows = TRUE, 
                                        rename_cols = TRUE) {
  
  tmp <- sample_map %>% select(sample_id, dap_dog_id)
  sample_map_sub <- tmp$dap_dog_id
  names(sample_map_sub) <- tmp$sample_id
  
  if (rename_rows) {
    if (sum(! feat_table$sample_id %in% names(sample_map_sub)) > 0) 
      stop("sample_id's in feature table do not match sample_map table")
    feat_table <- feat_table %>%
      mutate(sample_id = sample_map_sub[sample_id]) %>%
      rename(dap_dog_id = sample_id)
  }
  
  if (rename_cols) {
    tmp <- sample_map_sub[names(sample_map_sub) %in% colnames(feat_table)]
    feat_table <- feat_table %>% rename_at(vars(names(tmp)), ~tmp)
  }
  
  return(feat_table)
}






# Convert a Kraken table to convenient feature table. Effectively transposes a table.
reorganize_kraken_table <- function(df) {
  require(dplyr)
  require(tibble)
  
  name_of_first_col <- names(df)[1]
  
  df %>%
    column_to_rownames(name_of_first_col) %>%
    t() %>%
    data.frame() %>%
    rownames_to_column('sample_id')
}

# Read tab-delimited files
read_tsv2 <- function(p, cmnt = "") { 
  return(read_delim(p, comment = cmnt, delim = "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)) 
}

# Normalizes values to relative abundances.
# Assumes first column is sample id's and all the rest are taxa abundances.
convert_to_relab <- function(df) {
  require(vegan)
  relab_df <- decostand(df %>% select(-1), method = 'total', MARGIN = 1)
  relab_df <- bind_cols(df %>% select(1), relab_df)
  return(relab_df)
}

#' Remove rare features (columns) based on either a prevalence cutoff (i.e. % of
#'    non-zero values) or a mean-abundance cutoff, or both (in which case they 
#'    will be applied consecutively). 
#'    
#' @param feat_table Any feature table representing abundances (taxonomic/
#'    functional...)
#' @param filter_method One of 'prevalence', 'avg_abundance', 'both'
#' @param prevalence_cutoff The minimum portion of samples that should have a 
#'    non-zero value, to qualify the feature as non-rare (between 0 and 1, 
#'    default: 0.1). 
#' @param avg_abundance_cutoff Features with an average abundance (over samples)
#'    less than this threshold are considered rare and discarded.
#' 
#' @return An updated feature table 
remove_rare_features <- function(feat_table, 
                                 filter_method = 'both', 
                                 prevalence_cutoff = 0.1, 
                                 avg_abundance_cutoff = 0.005) {
  # Required libraries loading
  require(dplyr)
  
  # Argument verifications 
  if (prevalence_cutoff < 0 | prevalence_cutoff > 1) 
    error('Provided an invalid prevalence_cutoff value')
  if (avg_abundance_cutoff < 0) 
    error('Provided an invalid avg_abundance_cutoff value')
  
  # Initialize table with all current features
  new_feat_table <- feat_table
  n_samples <- nrow(feat_table)
  n_taxa_before_filter <- ncol(feat_table)-1
  
  # Prevalence calculations (number of non-zero values per feature)
  if (filter_method %in% c('prevalence', 'both')) {
    frequencies <- colSums(new_feat_table[,-1]>0) / n_samples
    new_feat_table <- new_feat_table[,c(TRUE, frequencies > prevalence_cutoff)]
  }
  
  # Average abundance calculations
  if (filter_method %in% c('avg_abundance', 'both')) {
    avg_abundances <- colSums(new_feat_table[,-1]) / n_samples
    new_feat_table <- new_feat_table[,c(TRUE, avg_abundances > avg_abundance_cutoff)]
  }
  
  n_taxa_after_filter <- ncol(new_feat_table)-1
  message(n_taxa_before_filter-n_taxa_after_filter, ' of ',n_taxa_before_filter, ' features were removed.')
  
  return(new_feat_table)
}

#' Run PERMANOVA (significance test for group-level differences) with adonis2 
#' (vegan library). 
#' 
#'    
#' @param disct_metadata A table with the metadata features to analyze. A 
#'    'sample id' column is expected, columns (metadata_cols_to_analyze) should 
#'    be categorical (or NA).
#' @param dist_matrix Microbiome distance matrix (column names correspond to 
#'    sample id's in the first column and are in the same order).
#' @param sample_id_column The name of the sample_id column. Present in both 
#'    `dist_matrix` and `disct_metadata`
#' @param metadata_cols_to_analyze The discrete columns in the disct_metadata 
#'    to run permanova.
#' @param metadata_cols_fix_effects A vector of columns to add as fixed 
#'   effects (e.g. age, sex etc.). 
#' @param fdr_threshold FDR threshold used to define significant findings. PCoA 
#'    plots will only be generated for significant findings.
#' @param plot_labels PCoA plot labels in case the discrete metadata are numeric values.
#' @param pcoa_fill optional fill color to the PCoA (should be as the same size
#'   of n_distinct(disct_metadata$metadata_cols_to_analyze)). If not provided, 
#'   the PCoA colors will be the defaults.
#' @param ord New order to the x-axis in the PCoA plot (should be the same values
#'    as in the disct_metadata$metadata_col_to_analyze).
#' @param add_centroid TRUE (default) if you wand to enter ellipse and centroid center
#'    to PCoA plot.
#' @param level The level at which to draw an ellipse (default 0.9)
#' @param pairs_mode if TRUE will generate multiple PCoA plots of each pair of categories
#'    in each column in each `metadata_cols_to_analyze`. It might be usefull tool for 
#'    additional visualization of the permanova results (the permanova will still run 
#'    on all the data). Default: FALSE (color all the data).
#' @return A list of three elements:
#'    `results` Table of p-values and FDRs for all `metadata_cols_to_analyze`.
#'    `plots` is a list of relevant plots: PCoA colored by discrete feature for all
#'    significant results.
#'    `pcoa_values` data frame of the pcoa values and the metadata.
get_permanova <- function(disct_metadata, 
                          dist_matrix, 
                          dist_metric_name,
                          sample_id_column = 'sample_id', 
                          metadata_cols_to_analyze = NULL, 
                          metadata_cols_fix_effects = NULL,
                          fdr_threshold = 0.1,
                          permutations = 9999,
                          plot_labels = NULL,
                          pcoa_fill = NULL,
                          ord = NULL,
                          add_centroid = TRUE,
                          seed = 111,
                          level = 0.9,
                          pairs_mode = FALSE
){
  
  # Required libraries loading
  require(dplyr)
  require(tibble)
  require(vegan)
  
  set.seed(seed)
  
  # Define metadata_cols_to_analyze as all columns if not provided, and check if exist if provided.
  if (!is.null(metadata_cols_to_analyze)) {
    if (!all(metadata_cols_to_analyze %in% colnames(disct_metadata)))
      stop("One of metadata_cols_to_analyze not found in metadata")
  } else {
    metadata_cols_to_analyze <- colnames(disct_metadata) 
    metadata_cols_to_analyze <- metadata_cols_to_analyze[!metadata_cols_to_analyze %in% c(sample_id_column)]
    warning("metadata_cols_to_analyze not provided, all metadata columns will be taken except sample_id_column.\n")
  }
  
  # More data verifications
  if (!sample_id_column %in% colnames(dist_matrix))
    stop("sample_id_column not found in dist_matrix")
  if (!sample_id_column %in% colnames(disct_metadata))
    stop("sample_id_column not found in metadata")
  if (!all(metadata_cols_to_analyze %in% colnames(disct_metadata)))
    stop("One of metadata_cols_to_analyze not found in metadata")
  if (nrow(disct_metadata) != nrow(dist_matrix))
    warning("# of samples is different between the distance matrix and the metadata. Only samples with both features and metadata will be used for the following analysis.\n")
  
  # Make sure all random/fixed effects are indeed in data
  if (!is.null(metadata_cols_fix_effects)) {
    if (!all(metadata_cols_fix_effects %in% colnames(disct_metadata)))
      stop("One of metadata_cols_fix_effects not found in metadata")
  }
  
  # Merge the dataframes
  df <- merge(x = dist_matrix, y = disct_metadata, by = sample_id_column, all.x=FALSE)  %>%
    remove_rownames %>% 
    column_to_rownames(var=sample_id_column)
  
  if (nrow(df) == 0) stop('No common samples found.')
  
  # Add fixed effects to a string that later will be added to the formula and check that they are not one of the metadata_cols_to_analyze.
  mixed_string = ""
  if (!is.null(metadata_cols_fix_effects)) {
    for (fixed_col in metadata_cols_fix_effects) {
      # Validate that all metadata_cols_fix_effects not in metadata_cols_to_analyze.
      if (fixed_col %in% metadata_cols_to_analyze)
        stop(paste0("Shared column between metadata_cols_fix_effects and metadata_cols_to_analyze: ", 
                    fixed_col, "\nTry again with different metadata_cols_fix_effects!\n"))
      
      # Remove rows with missing fixed effect values
      df <- df %>% filter(!is.na(get({{fixed_col}})))
      mixed_string =  paste0(mixed_string, " + tmp_df[,'", fixed_col, "'] ")
    }
  }
  
  # TODO: add random effects to formula!
  
  # Run adonis2
  pvals <- c()
  for (col in metadata_cols_to_analyze) {
    message('Working on variable: ', col)
    tmp_df <- df[!is.na(df[,col]), ]
    
    # Check that `col` is not a constant
    if (n_distinct(tmp_df[[col]]) == 1) stop(paste('Column', col, 'has only one value, cannot run PERMANOVA'))
    
    formula_string <- paste("tmp_df[, rownames(tmp_df)] ~ tmp_df[[col]]", mixed_string)
    permanova <- adonis2(as.formula(formula_string),
                         permutations = permutations)
    pvals[[col]] <- permanova$`Pr(>F)`[1]
  }
  
  # Define pcoa_values varible as NULL, in order to return it when there is no significant results.
  pcoa_values <- NULL
  
  # Call plot_colored_pcoa_discrete for significant features.
  all_fdrs <- p.adjust(pvals, 'fdr')
  plots <- list()
  if(!all(all_fdrs > fdr_threshold)) {
    for (metadata_feature in names(all_fdrs[all_fdrs <= fdr_threshold])){
      
      # Generate PCoA plots of all the data
      if (!pairs_mode){
        p12 <- plot_colored_pcoa_discrete(dist_matrix, disct_metadata, sample_id_column, metadata_feature, dist_metric_name, 1, 2, plot_labels, pcoa_fill, ord, add_centroid, level)
        plots[[paste0(metadata_feature,"_pc1&pc2")]] <- p12$p
        pcoa_values <- p12$pcoa_values
        p13 <- plot_colored_pcoa_discrete(dist_matrix, disct_metadata, sample_id_column, metadata_feature, dist_metric_name, 1, 3, plot_labels, pcoa_fill, ord, add_centroid, level)$p
        plots[[paste0(metadata_feature,"_pc1&pc3")]] <- p13
        p23 <- plot_colored_pcoa_discrete(dist_matrix, disct_metadata, sample_id_column, metadata_feature, dist_metric_name, 2, 3, plot_labels, pcoa_fill, ord, add_centroid, level)$p
        plots[[paste0(metadata_feature,"_pc2&pc3")]] <- p23
        
        # Generate PCoA plots for all pairs of comparisons (although the permanova is for all the data),
      } else {
        disct_metadata[, "metadata_feature"] <- disct_metadata[, metadata_feature]
        pairs <- combn(unique(disct_metadata$metadata_feature), 2)
        for (i in 1:ncol(pairs)) {
          pair <- pairs[, i]
          p12 <- plot_colored_pcoa_discrete(dist_matrix, disct_metadata, sample_id_column, metadata_feature, dist_metric_name, 1, 2, plot_labels, pcoa_fill, ord, add_centroid, level, pair)
          plots[[paste0(metadata_feature, "_pc1&pc2_", pair[1], "_", pair[2])]] <- p12$p
          pcoa_values <- p12$pcoa_values
          p13 <- plot_colored_pcoa_discrete(dist_matrix, disct_metadata, sample_id_column, metadata_feature, dist_metric_name, 1, 3, plot_labels, pcoa_fill, ord, add_centroid, level, pair)$p
          plots[[paste0(metadata_feature, "_pc1&pc3_", pair[1], "_", pair[2])]] <- p13
          p23 <- plot_colored_pcoa_discrete(dist_matrix, disct_metadata, sample_id_column, metadata_feature, dist_metric_name, 2, 3, plot_labels, pcoa_fill, ord, add_centroid, level, pair)$p
          plots[[paste0(metadata_feature, "_pc2&pc3_", pair[1], "_", pair[2])]] <- p23
          
        }
      }
    }
  }
  
  # Return a list of a data frame of the p-values and the FDRs, and the PCoA plots.
  results = data.frame(unlist(pvals), unlist(all_fdrs)) %>%
    rename(pvals = 1, FDR = 2) %>%
    tibble::rownames_to_column("feature")
  
  return(list(results = results, plots = plots, pcoa_values = pcoa_values))
}

#' Get genetic distance matrix
#' 
#' @param mat Any genetic relatedness matrix
#' @return A tidy genetic distance matrix (1-genetic relatedness matrix) with 0 along the diagonal
get_genetic_distance_matrix <- function(mat) {
  for (i in 1:nrow(mat)) {
    for (j in 1:ncol(mat)) {
      # zero values along diagonal
      if (i == j) {
        mat[i, j] <- 0
      } else {
        # zero negative values
        if(mat[i, j] < 0 ){
          mat[i,j] <- 0
        }
        #convert from similarity to distance matrix
        mat[i, j] <- 1 - mat[i, j] 
      }
    }}
  mat <- mat%>%
    as_tibble(rownames = "dap_dog_id")
  
  return(mat)
}

#' Run a Mantel test on two distance matrices
#'
#' @param distancesX A microbiome distance matrix (column names correspond 
#'    to sample id's, and first column also lists sample id's).
#' @param distancesY As above - for the second omic
#' @param sample_id_column The name of the sample_id column
#' @param corr_method Either "spearman", "kendall" or "pearson". Mantel 
#'    supports all 3.
#' @param n_permutations Number of permutations for significance test
#'
#' @return A list including a (1) ggplot object (mantel statistic) and (2) a 
#'    result table with statistic value and p-value.
#' @export
run_mantel_test <- function(distancesX, 
                            distancesY, 
                            sample_id_column = "sample_id",
                            corr_method = "spearman",
                            n_permutations = 499,
                            seed = 1234,
                            quiet = FALSE) {
  
  require(tibble)
  require(vegan)
  
  # Data verification 
  if (!sample_id_column %in% colnames(distancesX))
    stop("sample_id_column not found in distancesX")
  if (!sample_id_column %in% colnames(distancesY))
    stop("sample_id_column not found in distancesY")
  
  distX_samples <- distancesX[sample_id_column] %>%
    pull(sample_id_column)
  distY_samples <- distancesY[sample_id_column] %>%
    pull(sample_id_column)
  shared_samples <- intersect(distX_samples, distY_samples)
  
  if (length(shared_samples) == 0) stop("No shared samples between matrices")
  
  # Filter matrices to only contain shared samples
  if (length(shared_samples) != length(distX_samples) | length(shared_samples) != length(distY_samples)) {
    distancesX <- distancesX%>%
      filter(!!sym(sample_id_column) %in% shared_samples) %>%
      select(all_of(c(sample_id_column, shared_samples)))
    distancesY <- distancesY%>%
      filter(!!sym(sample_id_column) %in% shared_samples) %>%
      select(all_of(c(sample_id_column, shared_samples )))
    if (!quiet) message(str_c(length(shared_samples), " shared samples will be used for analysis\n"))
  }
  
  # Convert to matrix format
  distX_matrix <- distancesX %>%
    column_to_rownames(sample_id_column)
  
  distY_matrix <- distancesY%>%
    column_to_rownames(sample_id_column)
  
  # Make sure identical sample ordering
  distX_matrix <- distX_matrix[shared_samples, shared_samples]
  distY_matrix <- distY_matrix[shared_samples, shared_samples]
  
  # Make sure no NA's in matrices
  if (sum(is.na(distX_matrix)) > 0) stop("Found NAs in distance matrix X. No NAs allowed.")
  if (sum(is.na(distY_matrix)) > 0) stop("Found NAs in distance matrix Y. No NAs allowed.")
  
  # Calculate Mantel correlations
  set.seed(seed)
  mantel_results <- vegan::mantel(distX_matrix, distY_matrix, method = corr_method, permutations = n_permutations)
  
  result_table <- data.frame(p_value = mantel_results$signif, 
                             mantel_corr = mantel_results$statistic,
                             n_shared_samples = length(shared_samples))
  
  # Plot mantel statistic against shuffled versions
  p <- ggplot(data.frame(corr_stat = mantel_results$perm), aes(x = corr_stat)) +
    geom_histogram(color = 'black', fill = 'grey80', alpha = 0.4, bins = round(n_permutations / 10)) +
    geom_vline(xintercept = mantel_results$statistic, color = 'cadetblue3', linewidth = 2, alpha = 0.8) +
    theme_classic() +
    scale_y_continuous(expand = expansion(mult = c(0, .05))) +
    xlab(paste0('Mantel statistic (',corr_method,')')) +
    ylab('Count') +
    ggtitle('Mantel statistic - shuffled vs. true') +
    theme(plot.title = element_text(hjust = 0.5)) 
  
  return(list(results = result_table, plot = p))
}


#' Run a procrustes test on two omics, using distance matrices.
#'   The function returns a procrustes plot visualizing the superimposition and
#'   a p-value as returned from a permutation-based test (implemented in the 
#'   'vegan' package).
#' Note that significant p-values do not always imply "pretty" plots...
#'
#' @param distancesX A microbiome distance matrix (column names correspond 
#'    to sample id's, and first column also lists sample id's).
#' @param distancesY As above - for the second
#' @param omic_name_X Name of the first omic (for plot)
#' @param omic_name_Y Name of the second omic (for plot)
#' @param dist_metric Name of the distance metric used
#' @param sample_id_column The name of the sample_id column
#' @param n_permutations Number of permutations for significance test
#' @param quiet Print messages?
#' @param pcoa_n_dimensions Number of PCoA dimensions to use for procrustes 
#' @param seed For reproducability
#'
#' @return A list including a (1) ggplot object and (2) a p-value.
#' @export
#' 
run_procrustes_test <- function(distancesX, 
                                distancesY, 
                                omic_name_X = "omic X",
                                omic_name_Y = "omic Y",
                                dist_metric = "--", 
                                sample_id_column = 'sample_id', 
                                n_permutations = 499, 
                                quiet = FALSE,
                                pcoa_n_dimensions = 10,
                                transparency_reflecting_residuals = FALSE,
                                seed = 1234) {
  
  # Required libraries loading
  require(vegan)
  
  # Data validations
  # TODO
  
  # Extract only common samples
  common_samples <- intersect(distancesX[[sample_id_column]], distancesY[[sample_id_column]])
  if(!quiet) message('Found ', length(common_samples), ' samples common to both omics')
  d1 <- distancesX %>% filter(!!as.symbol(sample_id_column) %in% common_samples) %>% column_to_rownames(var = sample_id_column)
  d2 <- distancesY %>% filter(!!as.symbol(sample_id_column) %in% common_samples) %>% column_to_rownames(var = sample_id_column)
  d1 <- d1[common_samples, common_samples] 
  d2 <- d2[common_samples, common_samples] 
  
  # Calculate PCoA
  pcoa1 <- cmdscale(d1, k = min(pcoa_n_dimensions, nrow(d1)-1), eig = TRUE)
  pcoa2 <- cmdscale(d2, k = min(pcoa_n_dimensions, nrow(d2)-1), eig = TRUE)
  # Check out: ordiplot(pcoa1, type = "text", display = 'sites')
  
  # Procrustes test (test -> symmetric, rotation -> not symmetric)
  set.seed(seed)
  proc_res <- protest(pcoa1, pcoa2, scale = T, scores = "sites", permutations = n_permutations)
  
  # Extract rotated data, for plotting (also use scale factor)
  rotations_df <- data.frame(
    sample_id = common_samples,
    d1_ax1 = proc_res$X[,1],
    d1_ax2 = proc_res$X[,2],
    d2_rotated_ax1 = proc_res$Yrot[,1] / proc_res$scale,
    d2_rotated_ax2 = proc_res$Yrot[,2] / proc_res$scale,
    residual = residuals(proc_res)
  )
  
  # Generate plot
  title <- paste0('Procrustes plot (P-value: ', proc_res$signif, ')')
  subtitle <- paste('Based on',dist_metric,'distances')
  omic_name_Y <- paste(omic_name_Y,'(superimposed)',sep='\n')
  omic_colors <- c('grey20', 'darkred'); names(omic_colors) <- c(omic_name_X, omic_name_Y)
  
  if (transparency_reflecting_residuals) {
    p <- ggplot(rotations_df) +
      geom_point(aes(x = d1_ax1, y = d1_ax2, color = omic_name_X, alpha = 1/residual), size = 2.2) +
      geom_segment(aes(x = d1_ax1, y = d1_ax2, xend = d2_rotated_ax1, yend = d2_rotated_ax2, alpha = 1/residual), color=omic_colors[2], arrow = arrow(length = unit(0.2,"cm"))) +
      geom_point(aes(x = d2_rotated_ax1, y = d2_rotated_ax2, color = omic_name_Y, alpha = 1/residual), size = 2.2) +
      scale_color_manual(values = omic_colors, breaks = names(omic_colors)) +
      theme_classic() +
      ggtitle(title) +
      labs(subtitle = subtitle) +
      xlab('PCoA Axis 1') +
      ylab('PCoA Axis 2') +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(plot.subtitle = element_text(hjust = 0.5)) +
      theme(legend.title = element_blank())
  } else {
    p <- ggplot(rotations_df) +
      geom_point(aes(x = d1_ax1, y = d1_ax2, color = omic_name_X), alpha = 0.8, size = 2.2) +
      geom_segment(aes(x = d1_ax1, y = d1_ax2, xend = d2_rotated_ax1, yend = d2_rotated_ax2), alpha = 0.2, color=omic_colors[2], arrow = arrow(length = unit(0.2,"cm"))) +
      geom_point(aes(x = d2_rotated_ax1, y = d2_rotated_ax2, color = omic_name_Y), alpha = 0.8, size = 2.2) +
      scale_color_manual(values = omic_colors, breaks = names(omic_colors)) +
      theme_classic() +
      ggtitle(title) +
      labs(subtitle = subtitle) +
      xlab('PCoA Axis 1') +
      ylab('PCoA Axis 2') +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(plot.subtitle = element_text(hjust = 0.5)) +
      theme(legend.title = element_blank())
  }
  
  # Generate a plot of actual residuals compared to residuals from null models
  null_residuals <- data.frame()
  proc_shuf <- list()
  for(i in 0:10) {
    pcoa2_shuffled <- pcoa2
    
    # Shuffle all but 1st iteration
    if (i>0) {
      pcoa2_shuffled$points <-  pcoa2_shuffled$points[sample(1:nrow(pcoa2_shuffled$points)),]
      rownames(pcoa2_shuffled$points) <- sample(rownames(pcoa2_shuffled$points))
    }
    
    # Run procrustes
    proc_shuf[[as.character(i)]] <- procrustes(pcoa1$points, pcoa2_shuffled$points, scale = T, scores = "sites", symmetric = T)
    
    # Record residuals for plotting
    null_residuals <- bind_rows(
      null_residuals,
      data.frame(residual = residuals(proc_shuf[[as.character(i)]]), label = ifelse(i>0,'Shuffled','True'), shuffle_id = i)
    )
  }
  
  # Residuals plot
  p2 <- ggplot(null_residuals, aes(x = as.character(shuffle_id), y = residual, fill = label)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(alpha = 0.05, height = 0, width = 0.2) + 
    theme_classic() +
    ggtitle('Residuals in true vs. shuffled models') +
    ylab('Point-wise residuals') +
    xlab('True/shuffled procrustes models') +
    scale_fill_manual(values = c('True' = 'cadetblue3', 'Shuffled' = 'grey80')) +
    theme(legend.title = element_blank()) +
    theme(axis.text.x = element_blank()) +
    theme(plot.title = element_text(hjust = 0.5)) 
  
  # Generate a plot of the procrustes statistic compared to null models 
  # Note: Procrustes statistic = sqrt(1 - sum-of-squares-of-residuals)
  # Interpretation tip: The further away the "true" line is from the null values, the better.
  p3 <- ggplot(data.frame(t_stat = proc_res$t), aes(x = t_stat)) +
    geom_histogram(color = 'black', fill = 'grey80', alpha = 0.4, bins = 30) +
    geom_vline(xintercept = proc_res$t0, color = 'cadetblue3', linewidth = 2, alpha = 0.8) +
    theme_classic() +
    scale_y_continuous(expand = expansion(mult = c(0, .05))) +
    xlab('Procrustes statistic (sqrt(1-SSE))') +
    ylab('Count') +
    ggtitle('Procrustes statistic - shuffled vs. true') +
    theme(plot.title = element_text(hjust = 0.5)) 
  
  results_df = data.frame(
    p_value = proc_res$signif, 
    var_explained = proc_res$t0 ^ 2, # t_0 is the protest statistic, equivalent to a correlation 
    n_shared_samples = length(common_samples))
  
  list(results = results_df, plot = p, plot_residuals = p2, plot_procrustes_stat = p3)
}


#' Given a list of distance matrices, one of which defined as the response, MRM
#'   fits a linear regression modelling the 'response' distances as a linear 
#'   function of all other 'predictor' matrices. 
#' MRM original publication: Lichstein, J. 2007. Multiple regression on distance 
#'   matrices: A multivariate spatial analysis tool. Plant Ecology 188: 117-131.
#' MRM criticism: 
#'   https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0175194
#'
#' @param distances_list A named list of all distance matrices
#' @param response_matrix_name The name of the matrix to use as the response 
#'   matrix (the "y")
#' @param sample_id_column The name of the column holding sample ID's in each 
#'   matrix. Should be the same for all matrices.
#' @param n_permutations Number of MRM permutations (equivalent to Mantel's 
#'   test permutations)
#' @param quiet Control verbosity
#'
#' @return A named list with main MRM outputs
run_MRM <- function(distances_list, 
                    response_matrix_name, 
                    sample_id_column, 
                    n_permutations = 499,
                    method = 'spearman',
                    quiet = FALSE) {
  require(ecodist)
  
  # Verifications
  if (! response_matrix_name %in% names(distances_list))
    stop(response_matrix_name, " is missing from the distances_list")
  
  # Get samples common to all distance matrices
  common_samples <- Reduce(intersect, 
                           lapply(distances_list, 
                                  function(dm, sample_id_column) dm[[sample_id_column]], 
                                  sample_id_column = sample_id_column)
  )
  common_samples <- unname(common_samples)
  n_samps <- length(common_samples)
  if (n_samps == 0) stop('No samples shared across matrices')
  if (!quiet) message('Found ', n_samps, ' samples common to all omics')
  
  # Reorganize dist matrices to have same samples in the same order
  distances_aligned <- lapply(
    distances_list, 
    function(dm, sample_id_column, common_samples) {
      dm <- dm %>% 
        filter(!!as.symbol(sample_id_column) %in% common_samples) %>% 
        column_to_rownames(var = sample_id_column)
      dm <- dm[common_samples, common_samples]
      return(as.matrix(dm))
    }, 
    sample_id_column = sample_id_column,
    common_samples = common_samples
  )
  
  # Create formula for MRM
  predictor_names <- names(distances_list)[names(distances_list) != response_matrix_name]
  mrm_formula <- str_c(
    "as.dist(distances_aligned[['", 
    response_matrix_name, 
    "']]) ~ ",
    paste(sapply(
      predictor_names, 
      function(x) str_c("as.dist(distances_aligned[['", x, "']])")), 
      collapse = " + ")
  )
  
  # Run MRM
  mrm_results <- MRM(
    as.formula(mrm_formula), 
    nperm = n_permutations, 
    mrank = (method == 'spearman') # If FALSE, pearson is used
  ) 
  
  # Extract results
  mrm_coefs <- mrm_results$coef %>% 
    data.frame() %>% 
    rownames_to_column('matrix') %>% 
    rename(coef = 2) %>% 
    mutate(matrix = gsub('as\\.dist\\(distances_aligned\\[\\[\\"(.*)\\"\\]\\]\\)',"\\1", matrix)) 
  
  R2 <- unname(mrm_results$r.squared['R2'])
  AdjR2 <- 1 - ((1 - R2) * (n_samps - 1) / (n_samps - length(predictor_names) - 1))
  R2_pval = unname(mrm_results$r.squared['pval'])
  
  return(list(R2 = R2, AdjR2 = AdjR2, R2_pval = R2_pval, n_shared_samples = length(common_samples), coegs = mrm_coefs))
}

#' Train machine learning regression models based on all features in the 
#'    feature table, to assess predictability of each continuous metadata 
#'    feature. Model evaluation is carried out using a Spearman correlation 
#'    between predicted and actual values. 
#'    
#' @param feat_table A feature table (taxa abundances, pathways, alpha-diversity 
#'    metrics, etc.) with rows as samples and columns as features. A 'sample id' 
#'    column is expected in order to merge the table with the metadata.
#' @param metadata A table with the metadata features to analyze. A 'sample id' 
#'    column is expected.
#' @param sample_id_column The name of the sample_id column. Present in both 
#'    `feat_table` and `metadata`
#' @param vars_to_analyze A vector of columns names from the metadata 
#'    table for which the analysis should be carried out. 
#' @param ml_method One of "glm" (generalized linear model), or "rf" (random 
#'    forest regressor)
#' @param fdr_threshold FDR threshold used to define significant findings 
#' @param corr_threshold Correlation threshold used to define significant findings 
#' @param corr_method Correlation method for testing predicted values vs. true 
#'    values. Either 'pearson' or 'spearman'
#' @param n_ml_repeats Number of repeats for repeated cross-validation
#' @param n_ml_folds Number of folds for repeated cross-validation
#' @param feat_importance For RF models, 'altmann' method can be computed. Use 
#'   'none' to avoid feature importance (and decrease run times). 
#'   Default: 'none'.
#' @param rand_seed Random seed for reproducible results
#' 
#' @return A list with one element: 
#'    `results` contains summary results: correlations between predicted and 
#'    actual values, summarized over cv-folds.
#' @export
get_predictability_regression <- function(feat_table, 
                                          metadata, 
                                          sample_id_column, 
                                          vars_to_analyze, 
                                          ml_method = 'glm',
                                          fdr_threshold = 0.1,
                                          corr_threshold = 0.3,
                                          corr_method = "pearson",
                                          n_ml_repeats = 10,
                                          n_ml_folds = 5,
                                          feat_importance = 'none',
                                          rand_seed = 1111,
                                          quiet = FALSE) {
  # Required libraries loading
  require(rsample)
  require(dplyr)
  require(glmnet)
  require(ranger)
  require(metap) 
  
  # Argument validations 
  if(! all(vars_to_analyze %in% names(metadata)))
    stop("Invalid 'vars_to_analyze' argument. Some of the variables provided are missing from the metadata table")
  if(is.null(vars_to_analyze))
    stop("Invalid 'vars_to_analyze' argument. Cannot be NULL")
  if(! ml_method %in% c('glm', 'rf'))
    stop("Invalid 'ml_method' argument. Only 'glm' and 'rf' supported.")
  if (! corr_method %in% c('pearson', 'spearman'))
    stop("Invalid 'corr_method' argument. Only 'pearson' and 'spearman' supported.")
  
  # Main calculations
  n_samples <- nrow(feat_table)
  
  # Generate folds for cross-validation
  set.seed(rand_seed)
  folds <- vfold_cv(
    data.frame(sample_num = 1:n_samples),
    v = n_ml_folds,
    repeats = n_ml_repeats
  )
  
  # Collect evaluation stats here
  cv_correlations <- data.frame()
  cv_feat_imp <- data.frame()
  
  # Fix names of features (causes errors in some models)
  names(feat_table) <- make.names(names(feat_table))
  
  # Iterate over variables
  for (measure in vars_to_analyze) {
    if (!quiet) message("Training model for variable: ", measure)
    
    # Prepare temporary table for training (features + label only)
    tmp <- feat_table %>% 
      left_join(
        metadata[,c(sample_id_column, measure)] %>% 
          rename(Label = 2), 
        by = sample_id_column) %>%
      select(-any_of(c(sample_id_column)))
    
    # Cross validation loop
    for (i in 1:nrow(folds)) {
      fold_id <- paste(folds$splits[[i]]$id, collapse = " ")
      
      # Get train/test samples of this fold
      train_samples <- folds$splits[[i]]$in_id
      train_data <- tmp[train_samples,] %>% filter(!is.na(Label))
      test_data <- tmp[-train_samples,] %>% filter(!is.na(Label))
      
      # Train logistic regression / RF
      if (ml_method == 'glm') {
        lasso_model <- cv.glmnet(
          x = train_data %>% select(-Label) %>% as.matrix(), 
          y = train_data$Label, 
          grouped = (n_samples > 50)) # grouped = FALSE relevant for very small sample sets only
      } else if (ml_method == 'rf') {
        rf_model <- ranger(Label ~ ., data = train_data, importance = ifelse(feat_importance=='altmann','permutation','none'), seed = rand_seed)
        if (feat_importance != 'none') {
          feat_imp <- importance_pvalues(rf_model, method="altmann", formula = Label ~ ., data = train_data, seed = rand_seed)
        }
      }
      
      # Make prediction on held out samples
      if (ml_method == 'glm') {
        preds <- predict(lasso_model, 
                         test_data %>% select(-Label) %>% as.matrix(), 
                         s = "lambda.min") %>% unname()
      } else if (ml_method == 'rf') {
        preds <- predict(rf_model, test_data %>% select(-Label))$predictions
      }
      
      # Calculate model performance for this fold
      tmp_eval <- data.frame(preds = preds, true_label = test_data$Label)   
      cv_correlations <- bind_rows(
        cv_correlations,
        data.frame(
          measure = measure,
          fold_id = fold_id,
          test_cor = cor(tmp_eval$preds, tmp_eval$true_label, method = corr_method),
          test_cor_p = cor.test(tmp_eval$preds, 
                                tmp_eval$true_label, 
                                alternative = 'greater', # Check positive association
                                method = corr_method,
                                exact = FALSE)$p.value
        )
      )
      
      if (feat_importance != 'none')
        cv_feat_imp <- bind_rows(
          cv_feat_imp,
          as.data.frame(feat_imp) %>% 
            tibble::rownames_to_column('feature') %>%
            mutate(measure = measure, fold_id = fold_id)
        )
    }
  }
  
  sumlog2 <- function(ps) {ps <- ps[!is.na(ps)]; return(sumlog(ps)$p)}
  
  # Summarise over folds
  summary_correlations <- cv_correlations %>%
    group_by(measure) %>%
    summarise(n_folds = n(),
              mean_cor = mean(test_cor, na.rm = TRUE),
              sd_cor = sd(test_cor, na.rm = TRUE),
              p = sumlog2(test_cor_p),
              .groups = "drop") %>%
    mutate(fdr = p.adjust(p)) %>%
    mutate(Significant = fdr < fdr_threshold & mean_cor > corr_threshold)
  
  # Summarize feature importance (if requested by user)
  summary_feat_importance <- data.frame()
  if (feat_importance != 'none')
    summary_feat_importance <- cv_feat_imp %>%
    group_by(measure, feature) %>%
    summarise(n_folds = n(),
              mean_feat_importance = mean(importance, na.rm = TRUE),
              sd_feat_importance = sd(importance, na.rm = TRUE),
              p = sumlog2(pvalue),
              .groups = "drop") %>%
    mutate(fdr = p.adjust(p))
  
  if (! any(summary_correlations$Significant)) {
    message('No significantly well-predicted metadata measures')
  } else {
    message('The following measures were well-predicted by the microbiome using a ', toupper(ml_method), ' model:')
    message(paste(summary_correlations$measure[summary_correlations$Significant], collapse = ', '))
  }
  
  return(list(results = summary_correlations, feat_importance = summary_feat_importance))
}


#' Get mixed effect model: metadata_col ~ feature + fixed + (1 | random) for all:
#'  metadata_col in `metadata_cols_to_analyze` (from the metadata),
#'  feature in `feat_table` (except `sample_id_column`),
#'  fixed columns in the fixed_cols (from the metadata),
#'  random columns in the random_cols (from the metadata).
#'
#' @param feat_table A feature table with rows as samples and columns as features. 
#'    A 'sample id' column is expected in order to merge the table with the metadata.
##' @param metadata A table with the metadata features to analyze. A 'sample id' 
#'    column is expected.
#' @param metadata_cols_to_analyze A vector of columns names from the metadata 
#'    table for which the analysis should be carried out. If NULL, all available 
#'    columns will be analyzed.
#' @param sample_id_column The name of the sample_id column.
#' @param fixed_cols columns in the metadata for the fixed effect in the mixed effect model.
#' @param random_cols columns in the metadata for the random effect in the mixed effect model.
#' @param fdr_threshold FDR threshold used to define significant findings. Plots
#'    will only be generated for significant findings.
#' @param color_by_var column in the metadata to color the output scatter plots 
#'    by. By default resorts to the metadata columns being analyzed.
#'    
#' @return list of plots and results:
#'    plots- scatter plots of significant results, feature as x and metadata col in y.
#'    results- table of p-values and FDR for all the models.
#'
get_mixed_effect <- function(feat_table, 
                              metadata, 
                              metadata_cols_to_analyze = NULL, 
                              sample_id_column,
                              fixed_cols = NULL,
                              random_cols = NULL,
                              fdr_threshold = 0.1,
                              log_transform_y = FALSE,
                              color_by_var = NULL) {
  require(dplyr)
  require(tidyverse)
  require(lme4)
  require(ggpubr)
  require(ggplot2)
  require(ggpmisc)
  
  # Data verifications 
  if (!sample_id_column %in% colnames(feat_table))
    stop("sample_id_column not found in feat_tale")
  if (!sample_id_column %in% colnames(metadata))
    stop("sample_id_column not found in metadata")
  if (!is.null(metadata_cols_to_analyze)){
    if (!all(metadata_cols_to_analyze %in% colnames(metadata)))
      stop("One of metadata_cols_to_analyze not found in metadata")
  }
  if (nrow(metadata) != nrow(feat_table))
    warning("Note that nrow(metadata) != nrow(feat_table). Only samples with distance and metadata will be used for the following analysis.\n")
  if (!is.null(color_by_var))
    if (!color_by_var %in% colnames(metadata))
      stop("color_by_var not found in metadata")
  if (!is.null(fixed_cols)) {
    if (any(!fixed_cols %in% colnames(metadata)))
      stop("Some fixed_cols were not found in metadata")
  }
  if (!is.null(random_cols)) {
    if (any(!random_cols %in% colnames(metadata)))
      stop("Some random_cols were not found in metadata")
  }
  
  if (is.null(metadata_cols_to_analyze)){
    metadata_cols_to_analyze <- colnames(metadata) 
    metadata_cols_to_analyze <- metadata_cols_to_analyze[!metadata_cols_to_analyze %in% c(sample_id_column)]
    warning("metadata_cols_to_analyze not provided, all metadata columns will be taken except sample_id_column.")
  }
  
  # Convert vector of fixed effect and random effects to string.
  mixed_string = ""
  if (!is.null(fixed_cols)){
    for (fixed_col in fixed_cols)
      mixed_string =  paste(mixed_string, "+", fixed_col)
  }
  if (!is.null(random_cols)){
    for (random_col in random_cols)
      mixed_string =  paste(mixed_string,"+ ( 1 | ", random_col, ")")
  }
  
  # Merge the feat_table and the metadata data frames.
  df <- merge(x = feat_table, y = metadata, by = sample_id_column, all.x=FALSE)
  
  # Get the mixed effect model 
  results <- data.frame()
  for (metadata_col in metadata_cols_to_analyze){
    for (feature in colnames(feat_table[, -which(names(feat_table) == sample_id_column)])){
      
      # Organize table - remove rows with missing values
      df_tmp <- df %>% drop_na(metadata_col)
      df_tmp$metadata_col <- df_tmp[,metadata_col]
      if (log_transform_y) df_tmp$metadata_col <- log(df_tmp$metadata_col)
      df_tmp$feature <- df_tmp[,feature]
      
      # Create final LM formula
      formula_string <- paste("metadata_col", "~", "feature", mixed_string)
      
      if (is.null(random_cols)) {
        lmm1 <- lm(formula = as.formula(formula_string), data = df_tmp)
        coefs <- summary(lmm1)$coefficients
        print(coefs)
        pval <- ifelse('feature' %in% rownames(coefs), coefs['feature','Pr(>|t|)'], NA)
        estimate = coefs%>%as_tibble(rownames = "feature")%>%filter(feature == "feature")%>%pull(Estimate)
      } else {
        lmm1 <- lmer(formula = as.formula(formula_string), data = df_tmp)
        for_p_values <- parameters::p_value(lmm1)
        coefs <- summary(lmm1)$coefficients
        pval <- for_p_values$p[for_p_values$Parameter == "feature"] # ifelse('feature' %in% rownames(coefs), coefs['feature','Pr(>|z|)'], NA)
        estimate = coefs%>%as_tibble(rownames = "feature")%>%filter(feature == "feature")%>%pull(Estimate)
      }
      
      results <- bind_rows(
        results,
        data.frame(metadata_col = metadata_col,
                   feature = feature,
                   pval = pval,
                   N_samples = nrow(df_tmp),
                   estimate = estimate,
                   `formula` = formula_string)
      )
    }
  }
  
  results <- results %>%
    group_by(metadata_col) %>%
    mutate(FDR = p.adjust(pval, 'fdr')) %>%
    ungroup()
  
  # Create plots for significant results.
  
  if (is.null(color_by_var)){ # New
    color_by_var <- "tmp_col"
    df$tmp_col <- 0
  }
  
  plots <- list()
  for (metadata_col in metadata_cols_to_analyze){
    for (feature in colnames(feat_table[, -which(names(feat_table) == sample_id_column)])){
      
      # Plot settings
      # color_by_var_ <- ifelse(is.null(color_by_var), metadata_col, color_by_var # OLD 
      feature_ <- feature
      if (grepl(";g__", feature, fixed=TRUE))
        feature_ <- paste("g__", strsplit(feature_, ";g__")[[1]][2], sep="")
      
      fdr <- results$FDR[results$feature == feature & results$metadata_col == metadata_col]
      
      if (!is.na(fdr) & fdr <= fdr_threshold) {
        p <- ggplot(df, aes(x = .data[[metadata_col]], y = .data[[feature]])) + 
          # TODO
          # geom_path(aes(group = .data[[color_by_var_]]),
          #           color = 'darkgrey',
          #           arrow = arrow(length = unit(4, "points"), type = "closed")) +
          geom_point(aes(fill = .data[[color_by_var]]), color = 'black', shape = 21, size = 2, alpha = 0.5) +  # was: color_by_var_
          theme_classic() +
          ylab(feature_) +
          theme(legend.position = "none") + 
          geom_smooth(method = 'lm', formula = y ~ x, se=FALSE, color="darkred") +
          labs(subtitle = paste('FDR:', round(fdr, 4))) 
        plots[[paste(feature, "_", metadata_col)]] <- p
      }
    }
  }
  
  return(list(plots = plots, all_results = results, sig_results = results %>% filter(FDR <= fdr_threshold) %>% arrange(FDR)))
}

#' Check for if there is significant change in the distance from a reference value
#' for specific column. The function check for each sample what is the mean distance 
#' (according to dist_matrix, BC or WUF) to all the sample in the reference group.
#' This function perform two tests:
#' Regression: fit a linear regression mean_distance ~ metadata_col_to_analyze + metadata_cols_fix_effects.
#'    The reference_group_in_col is categorical, so the linear model return
#'    multiple p-values. In order to return only one p-value, this module perform
#'    likelihood ratio test to check if the model with the metadata_col_to_analyze
#'    as fixed effect is more likely than the model without it.
#' Wilcoxon: perform wilcoxon test that check if the distribution of the mean distance
#'    from the reference group is difference between all the other groups.
#' The function return two dataframes of results according to these two test.
#' In addition the function return boxplot of the means with significance value
#' according to the wilcoxon test (after FDR!).     
#'
#' @param disct_metadata A table with the metadata features to analyze. A 
#'    'sample id' column is expected, columns (metadata_cols_to_analyze) should 
#'    be categorical (or NA).
#' @param dist_matrix Microbiome distance matrix (column names correspond to 
#'    sample id's in the first column and are in the same order).
#' @param dist_metric_name 
#' @param metadata_col_to_analyze The discrete columns in the disct_metadata to 
#'    check.
#' @param reference_group_in_col The value in `metadata_col_to_analyze` to compare
#'    to.
#' @param sample_id_column The name of the sample_id column. Present in both 
#'    `dist_matrix` and `disct_metadata`
#' @param metadata_cols_fix_effects A vector of columns to add as fixed 
#'   effects (e.g. age, sex etc.). 
#' @param fdr_threshold fdr threshold used to define significant findings.
#'     default 0.1.
#' @param colors optional fill color to the boxplot. Format:
#'     list("category" = "color", ... ). default NULL.
#' @param ord optional ordering to the x axis (list). default NULL.
#'
#' @return A list of two elements:
#'    `regression_results` Table of p-value from the regression
#'    `wilcoxon_results` Table of all the p-values and the FDRs from 
#'         the wilcoxon test.
#'    `plot` boxplot of mean distances of each category to the  
#'        `metadata_col_to_analyze`.
get_distance_from_reference_group <- function(disct_metadata, 
                                               dist_matrix, 
                                               dist_metric_name,
                                               metadata_col_to_analyze, 
                                               reference_group_in_col,
                                               sample_id_column = 'sample_id', 
                                               metadata_cols_fix_effects = NULL,
                                               metadata_cols_random_effects = NULL,
                                               fdr_threshold = 0.1,
                                               colors = NULL,
                                               ord = NULL){
  
  # Required libraries loading
  require(lmtest)
  require(reshape2)
  require(tidyr)
  require(dplyr)
  require(lme4)
  
  # Validations
  if (!sample_id_column %in% colnames(dist_matrix))
    stop("sample_id_column not found in dist_matrix")
  if (!sample_id_column %in% colnames(disct_metadata))
    stop("sample_id_column not found in metadata")
  if (!(metadata_col_to_analyze %in% colnames(disct_metadata)))
    stop("metadata_col_to_analyze not found in metadata")
  if (nrow(disct_metadata) != nrow(dist_matrix))
    warning("# of samples is different between the distance matrix and the metadata. Only samples with both features and metadata will be used for the following analysis.\n")
  if(!any(grepl(reference_group_in_col, disct_metadata[[metadata_col_to_analyze]])))
    stop("reference_group_in_col not in metadata_col_to_analyze.")
  
  d <- dist_matrix %>%
    column_to_rownames(var=sample_id_column)
  d$row <- rownames(d)
  
  # Wide to long 
  long <- d %>% 
    pivot_longer(
      col = -row,
      names_to = "col",
      values_to = "distance"
    )
  
  disct_metadata$row <- disct_metadata$sample_id
  
  # Add metadata_col_to_analyze data to the first and the second samples to keep reference_group_in_col at the row and delete from the col
  long <- long %>% 
    left_join(disct_metadata %>% 
                select(c(!!sample_id_column, !!metadata_col_to_analyze)), by = c("row" = sample_id_column)) %>%
    rename(metadata_col_to_analyze_row = !!metadata_col_to_analyze) %>%
    left_join(disct_metadata %>% 
                select(c(!!sample_id_column, !!metadata_col_to_analyze)), by = c("col" = sample_id_column)) %>%
    rename(metadata_col_to_analyze_col = !!metadata_col_to_analyze) %>%
    filter(metadata_col_to_analyze_row == reference_group_in_col & metadata_col_to_analyze_col != reference_group_in_col)
  
  # Get means 
  long <- long %>% 
    group_by(col) %>%
    summarize(mean_distance = mean(distance)) %>%
    left_join(disct_metadata, by = c("col" = sample_id_column))
  
  
  # Add random effect
  mixed_string <- ""
  if (!is.null(metadata_cols_random_effects)){
    for (random_col in metadata_cols_random_effects)
      mixed_string =  paste(mixed_string,"+ ( 1 | ", random_col, ")")
  }
  
  # Liner regression
  
  #fit full and reduces models
  formula_string_full <- paste("mean_distance", "~", metadata_col_to_analyze, " + ", paste(metadata_cols_fix_effects, collapse = " + "), mixed_string)
  formula_string_reduced <- paste("mean_distance", "~", paste(metadata_cols_fix_effects, collapse = " + "), mixed_string)
  if (is.null(metadata_cols_random_effects)) {
    model_full <- lm(formula = as.formula(formula_string_full), data = long)
    model_reduced <- lm(formula = as.formula(formula_string_reduced), data = long)
  } else {
    model_full <- lmer(formula = as.formula(formula_string_full), data = long)
    model_reduced <- lmer(formula = as.formula(formula_string_reduced), data = long)
  }
  
  #perform likelihood ratio test for differences in models
  lr <- lrtest(model_full, model_reduced)
  pval <- lr$`Pr(>Chisq)`[2]
  
  # Order the x axis of the plot according to the following order
  long[, "metadata_col_to_analyze"] <- long[, metadata_col_to_analyze]
  if (!is.null(ord))
    long$metadata_col_to_analyze <- factor(long$metadata_col_to_analyze, levels = ord)
  
  # Add number of samples in each group to the x-axis ticks.
  xlabs <- list()
  for (val in as.character(unique(long$metadata_col_to_analyze)))
    xlabs[[val]] = paste0(val,"\n(", dim(long[as.character(long$metadata_col_to_analyze) == val,])[1], ")")
  
  # Plot boxplot of mean distance from the reference group
  p <- ggplot(data = long, 
              aes(x = metadata_col_to_analyze, y = mean_distance, fill = metadata_col_to_analyze)) +
    geom_boxplot() +
    geom_jitter(width = 0.2,  alpha = 0.3) +
    
    ylab(paste("Mean", "distance from", reference_group_in_col, "group")) +
    xlab(metadata_col_to_analyze) +
    # ggtitle(paste(dist_metric_name, ", log ratio test p-value:", plot_pval)) + 
    theme_classic() +
    theme(legend.position = "none")  +
    scale_x_discrete(labels=xlabs)
  
  # Add colors if specified
  if (!is.null(colors))
    p <- p + scale_fill_manual(values = colors)
  
  regression_df =  as.data.frame(list("metadata_col_to_analyze" = metadata_col_to_analyze,
                                      "reference_group_in_col" =  reference_group_in_col,
                                      "pvalue" = pval))
  
  # Perform wilcoxon test per each pair of unique values and add significance value
  uni_values <- distinct(long, !!sym(metadata_col_to_analyze), .keep_all= TRUE) %>%   pull({{ metadata_col_to_analyze }})
  wilcoxon_df = data.frame()
  all_pairs <- combn(uni_values, 2)
  for (i in 1:ncol(all_pairs)) {
    pair <- all_pairs[, i]
    pair_a <- long %>% filter(!!sym(metadata_col_to_analyze) == pair[1]) %>% pull(mean_distance)
    pair_b <- long %>% filter(!!sym(metadata_col_to_analyze) == pair[2]) %>% pull(mean_distance)
    wilcox_res <- wilcox.test(pair_a, pair_b, alternative = "two.sided")
    wilcoxon_df <- rbind(wilcoxon_df, data.frame(group1=c(pair[1]), group2=c(pair[2]),pval=c(wilcox_res$p.value)))
  }
  wilcoxon_df$FDR <- p.adjust(wilcoxon_df$pval, method = "fdr")
  wilcoxon_df <- wilcoxon_df %>% filter(FDR < fdr_threshold)  
  
  makeStars <- function(x){
    stars <- c("****", "***", "**", "*")
    vec <- c(0, 0.0001, 0.001, 0.01 ,0.05)
    i <- findInterval(x, vec)
    stars[i]
  }
  wilcoxon_df$star_fdr <- makeStars(wilcoxon_df$FDR)
  p <- p + stat_pvalue_manual(wilcoxon_df,
                              y.position = max(long$mean_distance),
                              step.increase = 0.1 ,
                              label = "star_fdr")
  
  return(list(regression_results =  regression_df, wilcoxon_resuls = wilcoxon_df, plot = p))
}


set.seed(123)

tmp_sample_map <- sample_map %>% filter(!is.na(weight_from_vet_or_survey_lbs)) %>%
  filter(life_stage_at_stool_collection != "Puppy") %>%
  mutate(new_life_stage_at_stool_collection = 
           ifelse(life_stage_at_stool_collection == "Young Adult" & 
                    runif(n()) <= 0.2, "Reference group young Adult", 
                  life_stage_at_stool_collection)
  )


run_continuous_maaslin <- function(feat_table, 
                                   metadata, 
                                   sample_id_column = 'sample_id', 
                                   metadata_cols_to_analyze = NULL,
                                   metadata_cols_fix_effects = NULL,
                                   metadata_cols_rand_effects = NULL,
                                   fdr_threshold = 0.1,
                                   output_dir = "MaAsLin2_output",
                                   return_maaslin_plots = FALSE,
                                   keep_maaslin_plots = FALSE,
                                   normalization_method = "NONE",
                                   transformation_method = "LOG",
                                   analysis_method = "LM",
                                   category_name_map = NULL,
                                   skip_fdr_correction = FALSE,
                                   logscale_y_for_plot = FALSE,
                                   quiet = FALSE) {
  # Required libraries loading
  require(Maaslin2)
  require(ggplot2)
  require(ggsignif)
  require(grid)
  require(png)
  require(dplyr)
  require(ggpubr)
  require(purrr)
  
  # Data verifications
  
  # In case maaslin plots required, change keep_maaslin_plots to TRUE.
  if (return_maaslin_plots & !(keep_maaslin_plots)){
    message("Can not return maaslin plots if keep_maaslin_plots = FALSE. keep_maaslin_plots was changed to TRUE.\n")
    keep_maaslin_plots = TRUE
  }
  
  # If no metadata_cols_to_analyze provided, run the MaAsLin2 on all the metadata columns.
  if (is.null(metadata_cols_to_analyze))
    metadata_cols_to_analyze = colnames(metadata)
  
  # Check if normalization method fits the data (if only one column: normalization method cannot be 'TSS')
  if ((ncol(feat_table) ==  2) & (normalization_method == "TSS")){
    stop("TSS normlaization is not allowed in case there is only one feature in feat_table. 
         Try again with normalization_method = 'NONE'.")
  }
  
  # Create output directory (will be deleted if not required).
  if (!dir.exists(output_dir))
    dir.create(output_dir)
  
  # Convert sample_id column to rownames for MaAsLin2 run.
  metadata <- metadata %>% remove_rownames %>% column_to_rownames(var=sample_id_column)
  feat_table <- feat_table %>% remove_rownames %>% column_to_rownames(var=sample_id_column)
  
  # Remove samples not in both tables, and order samples identically
  samples_to_include <- intersect(row.names(metadata), row.names(feat_table))
  metadata <- metadata[samples_to_include,]
  feat_table <- feat_table[samples_to_include,,drop=FALSE]
  
  # TODO: warn if some samples were not in both tables
  
  # Check that all the metadata_cols_rand_effects are not in metadata_cols_to_analyze.
  for (rand_col in metadata_cols_rand_effects){
    if (rand_col %in% metadata_cols_to_analyze)
      stop(paste0("Shared column between metadata_cols_rand_effects and metadata_cols_to_analyze: ", 
                  rand_col, "\nTry again with different metadata_cols_rand_effects!\n"))
  }
  
  # Go over fixed effects and perform one-hot encoding if needed
  fixed_cols_final_names <- c()
  for (fixed_col in metadata_cols_fix_effects) {
    
    # Validate that metadata_cols_fix_effects are not in metadata_cols_to_analyze.
    if (fixed_col %in% metadata_cols_to_analyze)
      stop(paste0("Shared column between metadata_cols_fix_effects and metadata_cols_to_analyze: ", 
                  fixed_col, "\nTry again with different metadata_cols_fix_effects!\n"))
    
    # Check if fixed column is numeric with less than 5 unique values or categorical. Perform one-hot encoding in these cases.
    if (((is.numeric(metadata[, fixed_col])) & 
         (n_distinct(metadata[, fixed_col]) < 5)) | 
        (!is.numeric(metadata[, fixed_col]))) {
      
      # Print a message if there are more than 5 categories in the metadata_cols_fix_effects because it might affect the FDR correction.
      if (n_distinct(metadata[[fixed_col]]) > 5)
        message("Too many unique values in ", fixed_col, " column. It might affect the FDR correction and omit some significant results!")
      
      # Perform one-hot encoding for columns if they have more then 2 categories.
      if (n_distinct(metadata[[fixed_col]]) > 2){
        message("Performing one-hot encoding on ", fixed_col)
        fixed_cols_final_names <- c(fixed_cols_final_names, paste(fixed_col, unique(metadata[[fixed_col]]), sep = '_'))
        tmp_metadata_rownames <- rownames(metadata)
        metadata <- metadata %>% mutate(value = 1)  %>% spread(fixed_col, value,  fill = 0, sep = '_')
        rownames(metadata) <- tmp_metadata_rownames
        
      } else { # Categorical with 2 values or less. One-hot is not needed 
        fixed_cols_final_names <- c(fixed_cols_final_names, fixed_col)
      }
      
    } else { # Numeric case, one-hot is not needed  
      fixed_cols_final_names <- c(fixed_cols_final_names, fixed_col)
    }
  }
  
  # Rum MaAsLin2
  maaslin2_signif_results <- data.frame() # Initialize
  maaslin2_all_results <- data.frame() # Initialize
  
  for (col in metadata_cols_to_analyze) {
    if (!quiet) message('Working on column: ', col)
    
    # Keep only the two selected values (also drop missing values)
    tmp_metadata <- metadata[(!is.na(metadata[[col]])),]# & (metadata[[col]] %in% pair), ]
    
    # Get the final list of samples for which DA will be calculated
    samples_to_include <- intersect(row.names(tmp_metadata), row.names(feat_table))
    tmp_metadata <- tmp_metadata[samples_to_include,,drop=F]
    tmp_feat_table <- feat_table[samples_to_include,,drop=F]
    
    # Remove constant columns if exist
    tmp_metadata <- tmp_metadata[, apply(tmp_metadata, 2, n_distinct) > 1]
    
    # Fix column names if needed
    colnames(tmp_metadata) <- make.names(colnames(tmp_metadata))
    if (!is.null(metadata_cols_rand_effects))
      make.names(metadata_cols_rand_effects)
    
    invisible(capture.output(tmp_results <- Maaslin2(
      input_data = tmp_feat_table, 
      input_metadata = tmp_metadata, 
      output = paste(output_dir, "/", col, sep=""),
      fixed_effects = make.names(c(col, fixed_cols_final_names)),
      random_effects = metadata_cols_rand_effects,
      normalization = normalization_method,
      transform = transformation_method,
      analysis_method = analysis_method,
      plot_scatter = keep_maaslin_plots,
      max_significance = fdr_threshold,
      plot_heatmap = FALSE
    )))
    
    #tmp_results$results$pair_a <- as.character(pair[[1]])
    #tmp_results$results$pair_b <- as.character(pair[[2]])
    
    # Save significant results only
    maaslin2_signif_results <- bind_rows(
      maaslin2_signif_results,
      tmp_results$results %>% filter(qval <= fdr_threshold)
    )
    
    # Save all results
    maaslin2_all_results <- bind_rows(
      maaslin2_all_results,
      tmp_results$results
    )
    #    }
  }
  
  # New FDR 
  if (skip_fdr_correction){
    maaslin2_all_results$FDR <- maaslin2_all_results$pval
  } else {
    maaslin2_all_results$FDR <- p.adjust(maaslin2_all_results$pval, 'fdr')
  }
  maaslin2_signif_results <- maaslin2_all_results[maaslin2_all_results$FDR < fdr_threshold,]
  
  # Return significant results only for column in metadata_cols_to_analyze (and not in metadata_cols_fix_effects).
  maaslin2_signif_results <- maaslin2_signif_results[maaslin2_signif_results$metadata %in% metadata_cols_to_analyze, ]
  
  # Generate plots to plots list
  plots <- list()
  
  if (nrow(maaslin2_signif_results) > 0) {
    if (return_maaslin_plots) {
      png_list <- fs::dir_ls(path = output_dir, recurse = TRUE, type = "file", glob = "*.png")
      for (png_dir in png_list){
        print(png_dir)
        plots[[png_dir]] <- rasterGrob(readPNG(png_dir))
      }
    } else {
      
      # Fix feature names to match how maaslin outputs them
      feat_table2 <- feat_table
      names(feat_table2) <- make.names(names(feat_table2))
      
      for (col in unique(maaslin2_signif_results$metadata)){
        for (feat in unique(maaslin2_signif_results$feature)){
          
          unique_signif_results <- maaslin2_signif_results[(maaslin2_signif_results$metadata == col) & (maaslin2_signif_results$feature == feat), ]
          if (dim(unique_signif_results)[1] == 0)
            next
          
          tmp_metadata <- metadata %>% drop_na(col)
          samples_to_include <- intersect(row.names(tmp_metadata), row.names(feat_table2))
          tmp_metadata <- tmp_metadata[samples_to_include,,drop=FALSE]
          feat_table2 <- feat_table[samples_to_include,,drop=FALSE]
          
          tmp <- data.frame(
            col_vals = tmp_metadata[[col]],
            feat_vals = feat_table2[[feat]]
          )
          
          
          FDR <- unique_signif_results[1,"FDR"][[1]]
          # Convert p-value to scientific style if needed 
          if (FDR < fdr_threshold / 1000){
            plot_FDR <- format(FDR, scientific = TRUE, digits=4)
          } else {
            plot_FDR <- round(FDR, 4)
          }
          
          p <- ggscatter(tmp, x = "col_vals", y = "feat_vals", title="MaAsLin2", subtitle = paste("FDR: ", plot_FDR),
                         color = "black", shape = 19, size = 3, alpha = 0.3, 
                         add = "reg.line",
                         add.params = list(color = "navy", fill = "lightgray"),
                         conf.int = TRUE,
                         cor.coeff.args = list(method = "spearman",  label.sep = "\n"),
                         xlab = col , ylab = feat
          ) + theme(plot.title = element_text(hjust = 0.5)) + font("xlab", size = 10) + font("ylab", size = 10) +
            font("xy.text", size = 10, color = "gray", face = "bold") 
          if (logscale_y_for_plot)
            p <- p + yscale("log10", .format = TRUE)
          p <- p + theme(plot.subtitle=element_text(size=10, hjust=0.5, color="black"))
          
          # p <- ggplot(tmp, aes(x = col_vals, group = col_vals, y = feat_vals)) +
          #   geom_point(color = 'black', alpha = 0.8, shape = 21, size = 4) +
          #   geom_smooth(method = 'lm', formula = feat_vals ~ col_vals, se=FALSE, color="darkred") +
          #   theme_classic() +
          #   xlab(col) +
          #   ylab(feat) +
          #   labs(subtitle = 'MaAsLin2') +
          #   theme(plot.subtitle = element_text(hjust = 0.5)) 
          # 
          # # Log scale y axis if needed
          # if (logscale_y_for_plot) {
          #   p <- p + yscale("log10", .format = TRUE)
          # }
          plots[[paste(col, feat, collapse = ';')]] <- p
        }
      }
    } 
  }
  
  if (! keep_maaslin_plots){
    unlink(output_dir, recursive=TRUE)
  }
  
  
  return(list(sig_results = maaslin2_signif_results, all_results = maaslin2_all_results, plots = plots))
}

#' Run MaAsLin2 (differential abundance) in case there are more than 2 discrete 
#' values in the col in metadata_cols_to_analyze.
#' This function will get the differential abundance of each pair and plot all 
#' the results in one boxplot.
#' In case you want to do "one vs. all" comparison you should doe one-hote encoding
#' to the specific column and run the get_differential_abundance function.
#'    
#' @param feat_table Any table with continuous features.
#' @param disct_metadata A table with the metadata features to analyze. A 
#'   'sample id' column is expected, columns (metadata_cols_to_analyze) should 
#'   be discrete.
#' @param sample_id_column The name of the sample id column.
#' @param metadata_cols_to_analyze The discrete columns in the disct_metadata 
#'   to run MaAaLin2.
#' @param metadata_cols_fix_effects A vector of columns to add as fixed 
#'   effects (e.g. age, sex etc.). 
#'   If numeric with less than 5 or categorical, will perform one-hot encoding.
#'   NULL to ignore.
#' @param metadata_cols_rand_effects A vector of columns to add as random 
#'   effects (e.g. subject idnetifier in the case of multiple samples per 
#'   subject, or batches in case of a multi-batch dataset). NULL to ignore.
#' @param fdr_threshold FDR threshold used to define significant findings. 
#'   Plots will only be generated for significant findings.
#' @param output_dir the directory to keep MaAsLin2 output.
#' @param return_maaslin_plots Boolean. FALSE (default) will return our own 
#'   plots instead of MaAsLin2 plots.
#' @param keep_maaslin_plots Boolean. FALSE (default) will not generate 
#'   MaAsLin2's default output files.
#' @param normalization_method A normalization method supported by MaAsLin2. 
#'   Defaults to 'NONE'. Other options are: CLR, CSS, NONE, TMM
#' @param transformation_method A transformation method supported by MaAsLin2. 
#'   Defaults to 'LOG'. Other options are: LOGIT, AST, NONE
#' @param analysis_method The model MaAsLin2 will run. Defaults to 'LM'.
#'   Other options are: CPLM, NEGBIN, ZINB
#' @param category_name_map raw category values will be mapped to names 
#'   according to this named vector, for plotting purposes. Also defines the 
#'   order of categories. NULL to skip.
#' @param logscale_y_for_plot log scale the Y axis in the plots (default: FALSE).
#' @param ord vector of values corresponding to metadata column to analyse unique
#'   values that will be the order of the x axis valus.
#' 
#' @return A list of two elements:
#'   `sig_results` contains statistical test results. Data frame with rows that 
#'    represent significant results. 
#'   `all_results` contains statistical test results. Data frame with rows that 
#'    represent all the results (provide the option to perform additional FDR 
#'    correction if needed,)
#'   `plots` is a list of relevant plots: boxplots that represent the 
#'    differential abundance of feature regarding the discrete metadata feature.
#'    Plots generated for adjusted p-value < fdr_threshold.
#'    The names of the plots list are the unique values (strings of the rows in 
#'    the results dataframe.)
get_differential_abundance <- function(feat_table, 
                                       disct_metadata, 
                                       sample_id_column = 'sample_id', 
                                       metadata_cols_to_analyze = NULL,
                                       metadata_cols_fix_effects = NULL,
                                       metadata_cols_rand_effects = NULL,
                                       fdr_threshold = 0.1,
                                       output_dir = "MaAsLin2_output",
                                       return_maaslin_plots = FALSE,
                                       keep_maaslin_plots = FALSE,
                                       normalization_method = "NONE",
                                       transformation_method = "LOG",
                                       analysis_method = "LM",
                                       category_name_map = NULL,
                                       skip_fdr_correction = FALSE,
                                       logscale_y_for_plot = FALSE,
                                       quiet = FALSE,
                                       ord=NULL) {
  # Required libraries loading
  require(Maaslin2)
  require(ggplot2)
  require(ggsignif)
  require(grid)
  require(png)
  require(dplyr)
  require(ggpubr)
  require(purrr)
  
  # Data verifications
  
  # In case maaslin plots required, change keep_maaslin_plots to TRUE.
  if (return_maaslin_plots & !(keep_maaslin_plots)){
    message("Can not return maaslin plots if keep_maaslin_plots = FALSE. keep_maaslin_plots was changed to TRUE.\n")
    keep_maaslin_plots = TRUE
  }
  
  # If no metadata_cols_to_analyze provided, run the MaAsLin2 on all the metadata columns.
  if (is.null(metadata_cols_to_analyze))
    metadata_cols_to_analyze = colnames(disct_metadata)
  
  # Check if metadata columns are discrete (1, 0 and NA).
  for (col in metadata_cols_to_analyze){ 
    if (length(unique(disct_metadata[,col])) > 3)
      warning(paste0("Too many values in disct_metadata in the ", col, " column"))
  }
  
  # Check if normalization method fits the data (if only one column: normalization method cannot be 'TSS')
  if ((ncol(feat_table) ==  2) & (normalization_method == "TSS")){
    stop("TSS normlaization is not allowed in case there is only one feature in feat_table. 
         Try again with normalization_method = 'NONE'.")
  }
  
  # Create output directory (will be deleted if not required).
  if (!dir.exists(output_dir))
    dir.create(output_dir)
  
  # Convert sample_id column to rownames for MaAsLin2 run.
  disct_metadata <- disct_metadata %>% remove_rownames %>% column_to_rownames(var=sample_id_column)
  feat_table <- feat_table %>% remove_rownames %>% column_to_rownames(var=sample_id_column)
  
  # Remove samples not in both tables, and order samples identically
  samples_to_include <- intersect(row.names(disct_metadata), row.names(feat_table))
  disct_metadata <- disct_metadata[samples_to_include,]
  feat_table <- feat_table[samples_to_include,,drop=FALSE]
  
  # TODO: warn if some samples were not in both tables
  
  # Check that all the metadata_cols_rand_effects are not in metadata_cols_to_analyze.
  for (rand_col in metadata_cols_rand_effects){
    if (rand_col %in% metadata_cols_to_analyze)
      stop(paste0("Shared column between metadata_cols_rand_effects and metadata_cols_to_analyze: ", 
                  rand_col, "\nTry again with different metadata_cols_rand_effects!\n"))
  }
  
  # Go over fixed effects and perform one-hot encoding if needed
  fixed_cols_final_names <- c()
  for (fixed_col in metadata_cols_fix_effects) {
    
    # Validate that metadata_cols_fix_effects are not in metadata_cols_to_analyze.
    if (fixed_col %in% metadata_cols_to_analyze)
      stop(paste0("Shared column between metadata_cols_fix_effects and metadata_cols_to_analyze: ", 
                  fixed_col, "\nTry again with different metadata_cols_fix_effects!\n"))
    
    # Check if fixed column is numeric with less than 5 unique values or categorical. Perform one-hot encoding in these cases.
    if (((is.numeric(disct_metadata[, fixed_col])) & 
         (n_distinct(disct_metadata[, fixed_col]) < 5)) | 
        (!is.numeric(disct_metadata[, fixed_col]))) {
      
      # Print a message if there are more than 5 categories in the metadata_cols_fix_effects because it might affect the FDR correction.
      if (n_distinct(disct_metadata[[fixed_col]]) > 5)
        message("Too many unique values in ", fixed_col, " column. It might affect the FDR correction and omit some significant results!")
      
      # Perform one-hot encoding for columns if they have more then 2 categories.
      if (n_distinct(disct_metadata[[fixed_col]]) > 2){
        message("Performing one-hot encoding on ", fixed_col)
        fixed_cols_final_names <- c(fixed_cols_final_names, paste(fixed_col, unique(disct_metadata[[fixed_col]]), sep = '_'))
        tmp_metadata_rownames <- rownames(disct_metadata)
        disct_metadata <- disct_metadata %>% mutate(value = 1)  %>% spread(fixed_col, value,  fill = 0, sep = '_')
        rownames(disct_metadata) <- tmp_metadata_rownames
        
      } else { # Categorical with 2 values or less. One-hot is not needed 
        fixed_cols_final_names <- c(fixed_cols_final_names, fixed_col)
      }
      
    } else { # Numeric case, one-hot is not needed  
      fixed_cols_final_names <- c(fixed_cols_final_names, fixed_col)
    }
  }
  
  # Rum MaAsLin2
  maaslin2_signif_results <- data.frame() # Initialize
  maaslin2_all_results <- data.frame() # Initialize
  
  for (col in metadata_cols_to_analyze) {
    if (!quiet) message('Working on column: ', col)
    
    # For multi-category variables, we run differential abundance separately for each pair of categories
    for (pair in combn(na.omit(unique(disct_metadata[[col]])), 2, simplify = FALSE)) {
      
      # Keep only the two selected values (also drop missing values)
      tmp_metadata <- disct_metadata[(!is.na(disct_metadata[[col]])) & (disct_metadata[[col]] %in% pair), ]
      
      # Get the final list of samples for which DA will be calculated
      samples_to_include <- intersect(row.names(tmp_metadata), row.names(feat_table))
      tmp_metadata <- tmp_metadata[samples_to_include,,drop=F]
      tmp_feat_table <- feat_table[samples_to_include,,drop=F]
      
      # Remove constant columns if exist
      tmp_metadata <- tmp_metadata[, apply(tmp_metadata, 2, n_distinct) > 1]
      
      # Fix column names if needed
      colnames(tmp_metadata) <- make.names(colnames(tmp_metadata))
      if (!is.null(metadata_cols_rand_effects))
        metadata_cols_rand_effects <- make.names(metadata_cols_rand_effects)
      
      invisible(capture.output(tmp_results <- Maaslin2(
        input_data = tmp_feat_table, 
        input_metadata = tmp_metadata, 
        output = paste(output_dir, "/", col, sep=""),
        fixed_effects = make.names(c(col, fixed_cols_final_names)),
        random_effects = metadata_cols_rand_effects,
        normalization = normalization_method,
        transform = transformation_method,
        analysis_method = analysis_method,
        plot_scatter = keep_maaslin_plots,
        max_significance = fdr_threshold,
        plot_heatmap = FALSE
      )))
      
      tmp_results$results$pair_a <- as.character(pair[[1]])
      tmp_results$results$pair_b <- as.character(pair[[2]])
      
      # Save significant results only
      maaslin2_signif_results <- bind_rows(
        maaslin2_signif_results,
        tmp_results$results %>% filter(qval <= fdr_threshold)
      )
      
      # Save all results
      maaslin2_all_results <- bind_rows(
        maaslin2_all_results,
        tmp_results$results
      )
    }
  }
  
  # New FDR 
  if (skip_fdr_correction){
    maaslin2_all_results$FDR <- maaslin2_all_results$pval
  } else {
    maaslin2_all_results$FDR <- p.adjust(maaslin2_all_results$pval, 'fdr')
  }
  maaslin2_signif_results <- maaslin2_all_results[maaslin2_all_results$FDR < fdr_threshold,]
  
  # Return significant results only for column in metadata_cols_to_analyze (and not in metadata_cols_fix_effects).
  maaslin2_signif_results <- maaslin2_signif_results[maaslin2_signif_results$metadata %in% metadata_cols_to_analyze, ]
  
  # Generate plots to plots list
  plots <- list()
  
  if (nrow(maaslin2_signif_results) > 0) {
    if (return_maaslin_plots) {
      png_list <- fs::dir_ls(path = output_dir, recurse = TRUE, type = "file", glob = "*.png")
      for (png_dir in png_list){
        print(png_dir)
        plots[[png_dir]] <- rasterGrob(readPNG(png_dir))
      }
    } else {
      
      # Fix feature names to match how maaslin outputs them
      feat_table2 <- feat_table
      names(feat_table2) <- make.names(names(feat_table2))
      
      
      for (col in unique(maaslin2_signif_results$metadata)){
        for (feat in unique(maaslin2_signif_results$feature)){
          
          unique_signif_results <- maaslin2_signif_results[(maaslin2_signif_results$metadata == col) & (maaslin2_signif_results$feature == feat), ]
          if (dim(unique_signif_results)[1] == 0)
            next
          
          
          tmp_metadata <- disct_metadata %>% tidyr::drop_na(col)
          samples_to_include <- intersect(row.names(tmp_metadata), row.names(feat_table2))
          tmp_metadata <- tmp_metadata[samples_to_include,,drop=FALSE]
          feat_table2 <- feat_table[samples_to_include,,drop=FALSE]
          
          tmp <- data.frame(
            col_vals = tmp_metadata[[col]],
            feat_vals = feat_table2[[feat]]
          )
          
          if (!is.null(category_name_map)) 
            tmp$col_vals <- factor(category_name_map[as.character(tmp$col_vals)],
                                   levels = unname(category_name_map))
          
          y_m <- max(tmp$feat_vals) 
          tmp$col_vals<- as.character(tmp$col_vals)
          
          # Add number of samples in each group to the x-axis ticks.
          xlabs <- list()
          for (val in as.character(unique(tmp$col_vals)))
            xlabs[[val]] = paste0(val,"\n(", dim(tmp[as.character(tmp$col_vals) == val,])[1], ")")
          
          # Sort the tmp dataset based on the specified order
          if (!is_null(ord)){     
            tmp$col_vals <- factor(tmp$col_vals, levels = ord)
            tmp <- tmp[order(tmp$col_vals), ]
          }
          
          p <- ggplot(tmp, aes(x = col_vals, group = col_vals, y = feat_vals)) +
            geom_boxplot(color = 'black', fill = 'lightblue', alpha = 0.7, outlier.shape = NA) +
            geom_jitter(height = 0, width = 0.1, alpha = 0.5, color = 'black', size = 2) +
            theme_classic() +
            xlab(col) +
            ylab(feat) +
            labs(subtitle = 'MaAsLin2') +
            theme(plot.subtitle = element_text(hjust = 0.5)) +
            scale_x_discrete(labels=xlabs)
          
          
          unique_signif_results$group1 <- unique_signif_results$pair_a
          unique_signif_results$group2 <- unique_signif_results$pair_b
          
          makeStars <- function(x){
            stars <- c("****", "***", "**", "*")
            vec <- c(0, fdr_threshold/1000, fdr_threshold/100, fdr_threshold/10 ,fdr_threshold)
            i <- findInterval(x, vec)
            stars[i]
          }
          
          unique_signif_results$FDR <- makeStars(unique_signif_results$FDR)
          
          p <- p + stat_pvalue_manual(unique_signif_results,
                                      y.position = max(tmp$feat_vals),
                                      step.increase = 0.1 ,
                                      label = "FDR")
          
          # Log scale y axis if needed
          if (logscale_y_for_plot) {
            p <- p + yscale("log10", .format = TRUE)
          }
          plots[[paste(col, feat, collapse = ';')]] <- p
        }
      }
    } 
  }
  
  if (! keep_maaslin_plots){
    unlink(output_dir, recursive=TRUE)
  }
  
  
  return(list(sig_results = maaslin2_signif_results, all_results = maaslin2_all_results, plots = plots))
}



#' Plot boxplot of values in feat_name (feat_table) grouped by 
#' metadata_col_to_analyze (disct_metadata).
#' Check if values feat_name in feat_table have different distribution in 
#' different groups of metadata_col_to_analyze in disct_metadata.
#' The different distribution is measured by a two-sided wilcoxon test. 
#' Significant annotation were were added according to fdr_threshold.
#' Significance annotations: *, <= fdr_threshold
#'                          **, <= fdr_threshold/10
#'                          ***, <= fdr_threshold/100
#'                          ****, <= fdr_threshold/100
#'
#' @param disct_metadata A table with the metadata features to analyze. A 
#'   'sample id' column is expected, columns (metadata_cols_to_analyze) should 
#'   be discrete.
#' @param feat_table A feature table (taxa abundances, pathways, alpha-diversity 
#'    metrics, etc.) with rows as samples and columns as features. A 'sample id' 
#'    column is expected in order to merge the table with the metadata.
#' @param feat_name The column in `feat_table` to analyse.
#' @param sample_id_column The name of the sample_id column. Present in both 
#'    `disct_metadata` and `feat_table`
#' @param metadata_col_to_analyze The metadata column to group by.
#' @param ord New order to the x-axis in the boxplot (should be the same values
#'    as in the disct_metadata$metadata_col_to_analyze).
#' @param colors Vector of colors to the boxplot (default: lightblue)
#' @param fdr_threshold FDR threshold used to define significant findings (used for
#'    the significant annotation above the boxplot).
#'
#' @return List of significant results and boxplots
feat_boxplot_wilcoxon <- function(disct_metadata, 
                                   feat_table, 
                                   feat_name,
                                   sample_id_column = 'sample_id', 
                                   metadata_col_to_analyze = NULL,
                                   ord = NULL,
                                   colors = "lightblue",
                                   fdr_threshold = 0.1){
  
  require("reshape2")
  require("dplyr")
  require("ggplot2")
  
  # Validations
  if(! metadata_col_to_analyze %in% names(disct_metadata))
    stop(metadata_col_to_analyze, " column is missing from metadata")
  if(! feat_name %in% names(feat_table))
    stop(feat_name, " column is missing from feat_table")
  if (!sample_id_column %in% colnames(feat_table))
    stop("sample_id_column not found in feat_table")
  if (!sample_id_column %in% colnames(disct_metadata))
    stop("sample_id_column not found in disct_metadata")
  if (nrow(disct_metadata) != nrow(feat_table))
    warning("# of samples is different between the distance matrix and the disct_metadata. Only samples with features and disct_metadata will be used for the following analysis.\n")
  
  # Keep only samples with values in metadata_col_to_analyze
  orig_nrows <- nrow(disct_metadata)
  disct_metadata <- disct_metadata %>% tidyr::drop_na(all_of(metadata_col_to_analyze))
  if (nrow(disct_metadata) < orig_nrows)
    message("Dropped ", orig_nrows - nrow(disct_metadata), " samples with missing ", metadata_col_to_analyze, " values.")
  
  df <- merge(x = feat_table, y = disct_metadata, by = sample_id_column, all.x=FALSE)
  df[,"metadata_col_to_analyze"] <- df[,metadata_col_to_analyze]
  df["feat_name"] <- df[, feat_name]
  
  # Perform thw wilcoxon test
  wilcoxon_df = data.frame()
  for (pair in combn(na.omit(unique(df$metadata_col_to_analyze)), 2, simplify = FALSE)){
    pair_a <- df[df$metadata_col_to_analyze == pair[1],]$feat_name
    pair_b <- df[df$metadata_col_to_analyze == pair[2],]$feat_name
    wilcox_res <- wilcox.test(pair_a, pair_b, alternative = "two.sided")
    wilcoxon_df <- rbind(wilcoxon_df, data.frame(a=c(pair[1]), b=c(pair[2]),pval=c(wilcox_res$p.value)))
  }
  
  # Perform FDR correction and keep only the significant values.
  wilcoxon_df$FDR <- p.adjust(wilcoxon_df$pval, 'fdr')
  wilcoxon_df <- wilcoxon_df[wilcoxon_df$FDR < fdr_threshold, ]
  
  # Order by default of `ord` if provided
  if (!is.null(ord))
    df$metadata_col_to_analyze <- factor(df$metadata_col_to_analyze , levels=ord)
  
  # Add number of samples in each group to the x-axis ticks.
  xlabs <- list()
  for (val in as.character(unique(df$metadata_col_to_analyze)))
    xlabs[[val]] = paste0(val,"\n(", dim(df[as.character(df$metadata_col_to_analyze) == val,])[1], ")")
  
  # Get the boxplots
  p <- ggplot(df, aes(x = metadata_col_to_analyze , 
                      group = metadata_col_to_analyze, 
                      y = feat_name)) +
    geom_boxplot(color = 'black', fill = colors, alpha = 0.7, outlier.shape = NA) +
    ylab(feat_name) +
    xlab(metadata_col_to_analyze) +
    theme_classic() + 
    labs(title=feat_name) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_discrete(labels=xlabs)
  
  # Add significance annotations.
  ym <- max(df[, feat_name])
  if (dim(wilcoxon_df)[1] > 0){
    
    mult_fact <- 0.1
    if (dim(wilcoxon_df)[1] > 4)
      mult_fact <- 0.05
    
    for (i in c(1:dim(wilcoxon_df)[1])){
      an = "*"
      if (wilcoxon_df$FDR[i] < (fdr_threshold / 10))
        an <- "**"
      if (wilcoxon_df$FDR[i] < (fdr_threshold / 100))
        an <- "***"
      if (wilcoxon_df$FDR[i] < (fdr_threshold / 1000))
        an <- "****"
      p <- p + ggsignif::geom_signif(annotations=an, y_position =ym+(mult_fact*i*ym), 
                                     xmin = wilcoxon_df$a[i], xmax =wilcoxon_df$b[i])
    }
  }
  return(list(results = wilcoxon_df, plot = p, df = df))
}


get_var_explained <- function(distance_mat, metadata, metadata_col, fixed_effects, random_effects = NULL) {
  
  # Get microbiome distance matrix and the metadata, to order them together.
  tmp <- merge(x = distance_mat, y = metadata, by = 'sample_id', all.x=FALSE)  %>%
    remove_rownames() %>% 
    column_to_rownames(var="sample_id")
  
  # Compute PCoA using the ape package
  dist_matrix <- as.dist(tmp[, rownames(tmp)])
  PCOA <- pcoa(dist_matrix)
  pcoa_values <- PCOA$vectors[,1:10] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("sample_id")
  
  # Add metadata to the pcoa values, will be used as the linear model input data frame.
  df <- pcoa_values %>% left_join(metadata, by = "sample_id")
  
  # Fit model (with PCs and confounders)
  regressors <- c(colnames(pcoa_values)[-1], fixed_effects)
  
  #Scale regressors
  numeric_regressors <- df%>%
    select_if(is.numeric) %>%
    colnames()%>%
    intersect(regressors)
  
  df <- df%>%
    mutate_at(vars(all_of(numeric_regressors)), datawizard::standardize)
  
  num_permutations <- 10
  random_permutations <- replicate(num_permutations, sample(regressors))
  
  final_pctexp <- data.frame("regressor" = c(sort(regressors),c("Residuals")))
  for (i in 1:num_permutations) {
    formula_string <- paste(metadata_col, "~", paste(random_permutations[, i], collapse = " + "))
    if(is.null(random_effects)){
      model <- lm(formula = as.formula(formula_string), data = df)
    }else{
      formula_string <- as.formula(str_c(formula_string, " + (1|", random_effects, ")"))
      model <- lmer(formula = as.formula(formula_string), data = df, REML = FALSE)
    }
    af <- anova(model)
    afss <- af$"Sum Sq"
    pctexp <- cbind(af,PctExp=afss/sum(afss)*100)
    pctexp <- pctexp %>% 
      select(PctExp) %>%
      mutate(regressor = row.names(.))
    colnames(pctexp)[colnames(pctexp) == "PctExp"] <- paste0("pctexp_perm", i)
    final_pctexp <- final_pctexp %>% left_join(pctexp, by = "regressor") 
  }
  
  # Remove the residuals from the final plot
  final_pctexp <- final_pctexp %>%
    filter(regressor != "Residuals")
  
  # Get mean and std for the final plot
  final_pctexp$mean_pctexp <- rowMeans(final_pctexp[,-1])
  final_pctexp$sd_pctexp <- apply(final_pctexp[,-c(1, (num_permutations+2))], 1, sd)
  
  # Order the x axis for the final plot
  final_pctexp$regressor  <- factor(final_pctexp$regressor  , levels=regressors)
  final_pctexp <- as_tibble(final_pctexp)
  
  # combine the microbiome PCs
  fixed_sum <- final_pctexp%>%
    filter(!str_detect(regressor, "Axis"))
  
  microbiome_sum <- final_pctexp%>%
    filter(str_detect(regressor, "Axis"))%>%
    summarize(across(where(is.numeric), sum))%>%
    mutate(regressor = "Microbiome")
  
  # orgenize data for ploting
  data_for_plot <- bind_rows(fixed_sum, microbiome_sum)%>%
    mutate(regressor = as.factor(regressor)%>%fct_reorder(mean_pctexp),
           metadata_col = metadata_col)
  
  return(data_for_plot)
}
