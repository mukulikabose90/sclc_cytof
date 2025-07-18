detect_batch_effect <- function (df, out_dir, downsample = NULL, norm_method = "scale", 
          xdim = 8, ydim = 8, seed = 382, markers = NULL, batch_col = "batch", 
          label_col = "label", name = "raw data") 
{
  cyCombine:::missing_package("outliers")
  cyCombine:::missing_package("Matrix")
  cyCombine:::missing_package("ggplot2", "CRAN")
  cyCombine:::missing_package("cowplot", "CRAN")
  cyCombine:::check_make_dir(out_dir)
  cyCombine:::check_colname(colnames(df), batch_col, location = "df")
  if (batch_col != "batch") {
    df$batch <- df[, batch_col]
  }
  df$batch <- as.factor(df$batch)
  if (length(levels(df$batch)) < 2) {
    stop("Error! Please provide a datasets with more than one batch.")
  }
  if (is.null(out_dir)) {
    stop("Error! Please speicify output directory.")
  }
  if (!is.null(downsample)) {
    message(paste0("Downsampling to ", downsample, " cells."))
    if (downsample <= nrow(df)) {
      set.seed(seed)
      df <- df %>% dplyr::slice_sample(n = downsample)
    }
    else {
      warning("Specified downsample parameter exceeds the number of cells in the dataset.")
    }
  }
  if (is.null(markers)) {
    markers <- df %>% cyCombine::get_markers()
  }
  if (!(label_col %in% colnames(df))) {
    message("Determining new cell type labels using SOM:")
    labels <- df %>% cyCombine::normalize(markers = markers, 
                                          norm_method = norm_method) %>% cyCombine::create_som(markers = markers, 
                                                                                               seed = seed, xdim = xdim, ydim = ydim)
    df$label <- labels
  }
  else {
    message("Using existig cell type labels.")
    df$label <- df[[label_col]]
  }
  if (length(levels(df$batch)) >= 3 & length(levels(df$batch)) <= 
      30) {
    emd <- df %>% dplyr::mutate(label = as.character(label)) %>% 
      cyCombine::compute_emd()
    markers_emd <- list()
    for (m in markers) {
      marker_emd <- lapply(emd, function(x) {
        x[[m]]
      })
      marker_emd <- lapply(marker_emd, function(x) {
        Matrix::forceSymmetric(x, uplo = "U")
      })
      marker_emd <- do.call(rbind, marker_emd)
      marker_emd <- Matrix::colMeans(as.matrix(marker_emd), 
                                     na.rm = T)
      markers_emd[[m]] <- marker_emd
    }
    marker_emd <- which(stats::p.adjust(sapply(markers_emd, 
                                               function(x) {
                                                 outliers::dixon.test(x, opposite = F)$p.value
                                               }), method = "BH") < 0.05 | stats::p.adjust(sapply(markers_emd, 
                                                                                                  function(x) {
                                                                                                    outliers::dixon.test(x, opposite = T)$p.value
                                                                                                  }), method = "BH") < 0.05)
    message(paste0("\nThere are ", length(marker_emd), " markers that appear to be outliers in a single batch:"))
    message(paste(markers[marker_emd], collapse = ", "))
    counts <- table(df$batch, df$label)
    perc <- (counts/Matrix::rowSums(counts)) * 100
    out_cl <- which(stats::p.adjust(apply(perc, 2, function(x) {
      outliers::dixon.test(x, opposite = F)$p.value
    }), method = "BH") < 0.05 | stats::p.adjust(apply(perc, 
                                                      2, function(x) {
                                                        outliers::dixon.test(x, opposite = T)$p.value
                                                      }), method = "BH") < 0.05)
    message(paste0("\nThere are ", length(out_cl), " clusters, in which a single cluster is strongly over- or underrepresented."))
    for (cl in out_cl) {
      message(paste0("The cluster percentages for each batch in cluster ", 
                     cl, " are:"))
      message(paste(names(perc[, cl]), "=", round(perc[, 
                                                       cl], 2), "%", collapse = ", "))
      exp_markers <- df %>% dplyr::filter(label == cl) %>% 
        dplyr::summarise(dplyr::across(cyCombine::get_markers(df), 
                                       stats::median)) %>% unlist() %>% sort(decreasing = T)
      message(paste0("The cluster expresses ", paste(names(exp_markers[which(exp_markers > 
                                                                               1)]), collapse = ", ")), ".")
      cat("\n")
    }
  }
  message("Making UMAP plots for up to 50,000 cells.")
  if (nrow(df) > 50000) {
    set.seed(seed)
    df <- df %>% dplyr::slice_sample(n = 50000)
  }
  plot_dimred_full(df, name, type = "umap", markers = NULL, 
                              seed = 473, out_dir)
  message(paste0("Saved UMAP plot for batches and labels here: ", 
                 out_dir, " as UMAP_batches_labels.png."))
  message(paste0("Saved UMAP plot colored by each marker in directory: ", 
                 out_dir, "/UMAP_markers.\n"))
  message("Done!")
}


plot_dimred_full <- function (df, name, type = "umap", markers = NULL, seed = 473, 
          out_dir = NULL) 
{
  if (type == "umap") 
    cyCombine:::missing_package("uwot", "CRAN")
  cyCombine:::missing_package("ggridges", "CRAN")
  cyCombine:::missing_package("grDevices", "CRAN")
  cyCombine:::missing_package("RColorBrewer", "CRAN")
  cyCombine:::missing_package("ggplot2", "CRAN")
  if (is.null(out_dir)) {
    stop("Error! Please speicify output directory.")
  }
  else {
    cyCombine:::check_make_dir(out_dir)
  }
  if (is.null(markers)) {
    markers <- cyCombine::get_markers(df)
  }
  Batch <- df$batch %>% as.factor()
  Label <- df$label %>% as.factor()
  df <- df %>% dplyr::select(dplyr::all_of(markers))
  if (type == "pca") {
    pca <- df %>% stats::prcomp(scale. = TRUE, center = TRUE)
    df <- cbind.data.frame(pca$x, as.factor(Batch), as.factor(Label), 
                           df)
    colnames(df)[3:ncol(df)] <- c("Batch", "Label", markers)
  }
  else if (type == "umap") {
    set.seed(seed)
    umap <- df %>% uwot::umap(n_neighbors = 15, min_dist = 0.2, 
                              metric = "euclidean")
    df <- cbind.data.frame(umap, as.factor(Batch), as.factor(Label), 
                           df)
    colnames(df) <- c("UMAP1", "UMAP2", "Batch", "Label", 
                      markers)
  }
  batch_plot <- df %>% ggplot2::ggplot(ggplot2::aes_string(x = colnames(df)[1], 
                                                           y = colnames(df)[2])) + ggplot2::geom_point(ggplot2::aes_string(color = "Batch"), 
                                                                                                       alpha = 0.3, size = 0.4, shape = 1) + ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(alpha = 1, 
                                                                                                                                                                                                               size = 1))) + ggplot2::theme_bw() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + 
    ggplot2::ggtitle(paste(toupper(type), "-", name))
  label_plot <- df %>% ggplot2::ggplot(ggplot2::aes_string(x = colnames(df)[1], 
                                                           y = colnames(df)[2])) + ggplot2::geom_point(ggplot2::aes_string(color = "Label"), 
                                                                                                       alpha = 0.3, size = 0.4, shape = 1) + ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(alpha = 1, 
                                                                                                                                                                                                               size = 1))) + ggplot2::theme_bw() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + 
    ggplot2::ggtitle(paste(toupper(type), "-", name))
  cyCombine::plot_save_two(batch_plot, label_plot, paste0(out_dir, 
                                                          "/UMAP_batches_labels.png"))
  cyCombine:::check_make_dir(paste0(out_dir, "/UMAP_markers"))
  
  marker_plots <- list()
  for (m in markers) {
    cat(m,"\n")
    p <- ggplot2::ggplot(df, ggplot2::aes_string(x = colnames(df)[1], y = colnames(df)[2])) + 
      ggplot2::geom_point(ggplot2::aes_string(color = m), alpha = 0.3, size = 0.4) + ggplot2::scale_color_gradientn(m, colors = (grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral"))))(50)) + 
      ggplot2::theme_bw() + 
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + 
      ggplot2::ggtitle(paste(toupper(type), "-", name))
    suppressMessages(ggplot2::ggsave(p, filename = paste0(out_dir, 
                                                          "/UMAP_markers/UMAP_batches_labels", m, ".png")))
  }
}


detect_batch_effect_express <- function (df, out_dir, batch_col = "batch", downsample = NULL, 
                                         seed = 472) 
{
  cyCombine:::missing_package("Matrix")
  cyCombine:::missing_package("ggridges", "CRAN")
  cyCombine:::missing_package("ggplot2", "CRAN")
  cyCombine:::missing_package("cowplot", "CRAN")
  cyCombine:::missing_package("grDevices", "CRAN")
  message("Starting the quick(er) detection of batch effects.")
  cyCombine:::check_make_dir(out_dir)
  cyCombine:::check_colname(colnames(df), batch_col, location = "df")
  if (batch_col != "batch") {
    df$batch <- df[, batch_col]
  }
  df$batch <- as.factor(df$batch)
  if (length(levels(df$batch)) < 2) {
    stop("Error! Please provide a datasets with more than one batch.")
  }
  df$label <- 1
  if (!("id" %in% names(df))) {
    df$id <- 1:nrow(df)
  }
  if (!is.null(downsample)) {
    message(paste0("Downsampling to ", downsample, " cells."))
    if (downsample <= nrow(df)) {
      set.seed(seed)
      df <- df %>% dplyr::slice_sample(n = downsample)
    }
    else {
      warning("Specified downsample parameter exceeds the number of cells in the dataset.")
    }
  }
  message("Making distribution plots for all markers in each batch.")
  all_markers <- df %>% cyCombine::get_markers()
  grDevices::pdf(NULL)
  p <- list()
  for (c in 1:length(all_markers)) {
    p[[c]] <- ggplot2::ggplot(df, ggplot2::aes_string(x = all_markers[c], 
                                                      y = "batch")) + ggridges::geom_density_ridges(ggplot2::aes(color = batch, 
                                                                                                                 fill = batch), alpha = 0.4, quantile_lines = TRUE) + 
      ggplot2::theme_bw()
  }
  suppressMessages(cowplot::save_plot(paste0(out_dir, "/distributions_per_batch.png"), 
                                      
                                      cowplot::plot_grid(plotlist = p, nrow = round(length(all_markers)/6)), 
                                      base_width = 30, base_height = length(all_markers)))
  message(paste0("Saved marker distribution plots here: ", 
                 out_dir, "/distributions_per_batch.png."))
  message("Applying global EMD-based batch effect detection.")
  emd <- df %>% dplyr::arrange(id) %>% cyCombine::compute_emd()
  emd_markers <- cbind.data.frame(sapply(emd[[1]], mean, na.rm = T), 
                                  sapply(emd[[1]], stats::sd, na.rm = T))
  emd <- lapply(emd[[1]], function(x) {
    Matrix::forceSymmetric(x, uplo = "U")
  })
  batch_means <- lapply(emd, function(x) {
    Matrix::colMeans(as.matrix(x), na.rm = T)
  })
  emd_markers <- cbind.data.frame(emd_markers, rownames(emd_markers))
  colnames(emd_markers) <- c("mean", "sd", "marker")
  emd_markers$marker <- factor(emd_markers$marker, levels = emd_markers$marker[order(emd_markers$mean, 
                                                                                     decreasing = T)])
  p <- ggplot2::ggplot(emd_markers, ggplot2::aes(x = marker, 
                                                 y = mean)) + ggplot2::geom_bar(stat = "identity", ggplot2::aes(fill = mean)) + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, 
                                                       vjust = 0.5, hjust = 1)) + ggplot2::ylab("Mean EMD") + 
    ggplot2::xlab("") + ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - 
                                                              sd, ymax = mean + sd), width = 0.2, position = ggplot2::position_dodge(0.9))
  suppressMessages(ggplot2::ggsave(p, filename = paste0(out_dir, 
                                                        "/emd_per_marker.png")))
  message(paste0("Saved EMD plot here: ", out_dir, "/emd_per_marker.png.\n"))
  any_outliers <- F
  for (m in emd_markers$marker[emd_markers$mean > stats::median(emd_markers$mean)]) {
    if (any(batch_means[[m]] > (stats::quantile(batch_means[[m]], 
                                                0.75) + stats::IQR(batch_means[[m]]) * 3))) {
      found_outliers <- which(batch_means[[m]] > (stats::quantile(batch_means[[m]], 
                                                                  0.75) + stats::IQR(batch_means[[m]]) * 3))
      message(paste0(m, " has clear outlier batch(es): ", 
                     paste(names(found_outliers), collapse = ", ")))
      summary_non <- df %>% dplyr::filter(!(batch %in% 
                                              names(found_outliers))) %>% dplyr::pull(m) %>% 
        summary()
      summary_out <- df %>% dplyr::filter(batch %in% names(found_outliers)) %>% 
        dplyr::pull(m) %>% summary()
      message("Summary of the distribution in the OUTLIER batch(es):")
      message(paste(names(summary_out), "=", round(summary_out, 
                                                   2), collapse = ", "))
      cat("\n")
      message("Summary of the distribution in the non-outlier batches:")
      message(paste(names(summary_non), "=", round(summary_non, 
                                                   2), collapse = ", "))
      cat("\n\n")
      any_outliers <- T
    }
  }
  if (!(any_outliers)) {
    message("None of the markers has very strong outlier batches, consult plots for more general deviations.")
  }
  message("Making MDS plot for median protein expression per sample.")
  median_expr <- df %>% dplyr::group_by(sample) %>% dplyr::summarise_at(cyCombine::get_markers(df), 
                                                                        stats::median)
  dist_mat <- as.matrix(stats::dist(median_expr[, all_markers]))
  rownames(dist_mat) <- colnames(dist_mat) <- median_expr$sample
  mds <- as.data.frame(stats::cmdscale(dist_mat, k = 2))
  mds$sample <- as.factor(rownames(mds))
  mds$batch <- as.factor(df$batch[match(mds$sample, df$sample)])
  p <- ggplot2::ggplot(mds, ggplot2::aes(V1, V2)) + ggplot2::geom_point(ggplot2::aes(colour = batch), 
                                                                        size = 2) + ggplot2::labs(x = "MDS1", y = "MDS2", title = "MDS plot") + 
    ggplot2::theme_bw()
  suppressMessages(ggplot2::ggsave(p, filename = paste0(out_dir, 
                                                        "/MDS_per_batch.png")))
  message(paste0("Saved MDS plot here: ", out_dir, "/MDS_per_batch.png"))
  message("Done!")
}
