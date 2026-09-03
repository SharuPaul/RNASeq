// Differential expression with DESeq2 

process merge_featurecounts {
    label 'deseq2'
    container 'quay.io/biocontainers/bioconductor-deseq2:1.50.2--r45ha27e39d_0'
    conda 'bioconda::bioconductor-deseq2=1.50.2'

    input:
    path(featurecounts_files)

    output:
    path('gene_count_matrix.tsv'), emit: count_matrix

    publishDir "${params.outdir}/06_DifferentialExpression", mode: 'copy'

    script:
    """
    Rscript - <<'RSCRIPT'
    input_files <- list.files(".", pattern = "_genecounts\\\\.txt\$", full.names = TRUE)

    if (length(input_files) == 0) {
      stop("No featureCounts gene count files matching '*_genecounts.txt' were found.")
    }

    read_featurecounts <- function(file) {
      tab <- read.delim(file, comment.char = "#", check.names = FALSE)

      if (!"Geneid" %in% names(tab)) {
        stop("Missing Geneid column in ", file)
      }

      sample <- sub("_genecounts\\\\.txt\$", "", basename(file))
      sample <- sub("_sorted\$", "", sample)
      raw_counts <- suppressWarnings(as.numeric(tab[[ncol(tab)]]))

      if (any(is.na(raw_counts))) {
        stop("Non-numeric count values found in ", file)
      }

      counts <- data.frame(
        Geneid = as.character(tab[["Geneid"]]),
        count = as.integer(round(raw_counts)),
        check.names = FALSE
      )
      names(counts)[2] <- sample
      counts
    }

    count_tables <- lapply(input_files, read_featurecounts)
    sample_names <- vapply(count_tables, function(tab) names(tab)[2], character(1))

    if (anyDuplicated(sample_names)) {
      duplicated_samples <- unique(sample_names[duplicated(sample_names)])
      stop("Duplicate sample names after suffix cleanup: ", paste(duplicated_samples, collapse = ", "))
    }

    merged_counts <- Reduce(
      function(left, right) merge(left, right, by = "Geneid", all = TRUE, sort = FALSE),
      count_tables
    )
    merged_counts[is.na(merged_counts)] <- 0

    write.table(
      merged_counts,
      file = "gene_count_matrix.tsv",
      sep = "\\t",
      quote = FALSE,
      row.names = FALSE
    )
    RSCRIPT
    """
}

process deseq2 {
    label 'deseq2'
    container 'quay.io/biocontainers/bioconductor-deseq2:1.50.2--r45ha27e39d_0'
    conda 'bioconda::bioconductor-deseq2=1.50.2'
    tag "$contrast"

    input:
    path(count_matrix)
    path(metadata)
    val(design_formula)
    val(contrast)

    output:
    path('deseq2_results_*.tsv'), emit: results
    path('deseq2_summary_*.tsv'), emit: summary
    path('deseq2_normalized_counts.tsv'), emit: normalized_counts
    path('deseq2_vst_counts.tsv'), emit: vst_counts
    path('deseq2_pca.pdf'), emit: pca
    path('deseq2_pca_ellipse.pdf'), emit: pca_ellipse
    path('deseq2_pca_centroid_distances.tsv'), emit: pca_distances
    path('deseq2_pca_centroids.tsv'), emit: pca_centroids
    path('deseq2_ma_plot_*.pdf'), emit: ma_plot
    path('deseq2_volcano_*.pdf'), emit: volcano
    path('deseq2_volcano_*.png'), emit: volcano_png
    path('deseq2_session_info.txt'), emit: session_info

    publishDir {
    def d = design_formula
        .replaceAll(/[^A-Za-z0-9]+/, '_')
        .replaceAll(/^_+|_+$/, '')

    "${params.outdir}/06_DifferentialExpression/${d}"
    }, mode: 'copy'

    script:
    """
    Rscript - "${count_matrix}" "${metadata}" "${design_formula}" "${contrast}" <<'RSCRIPT'
    suppressPackageStartupMessages({
      library(DESeq2)
      library(ggplot2)
    })

    args <- commandArgs(trailingOnly = TRUE)

    if (length(args) != 4) {
      stop("Expected arguments: count_matrix metadata design_formula contrast")
    }

    count_matrix_file <- args[[1]]
    metadata_file <- args[[2]]
    design_formula <- args[[3]]
    contrast_spec <- args[[4]]

    read_delimited <- function(file) {
      first_line <- readLines(file, n = 1, warn = FALSE)
      sep <- if (grepl("\\t", first_line, fixed = TRUE)) "\\t" else ","

      read.table(
        file,
        header = TRUE,
        sep = sep,
        quote = '"',
        comment.char = "",
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }

    contrast_entries <- trimws(strsplit(contrast_spec, ";", fixed = TRUE)[[1]])
    contrast_entries <- contrast_entries[contrast_entries != ""]

    if (length(contrast_entries) == 0) {
      stop("At least one contrast is required.")
    }

    parse_contrast <- function(contrast_entry) {
      contrast_parts <- trimws(strsplit(contrast_entry, ",", fixed = TRUE)[[1]])

      if (length(contrast_parts) != 3 || any(contrast_parts == "")) {
        stop("Each contrast must use the format 'variable,numerator,denominator'. Invalid contrast: ", contrast_entry)
      }

      contrast_parts
    }

    contrasts <- lapply(contrast_entries, parse_contrast)

    make_safe_name <- function(values) {
      safe_name <- paste(values, collapse = "_")
      safe_name <- gsub("[^A-Za-z0-9_.-]+", "_", safe_name)
      safe_name <- gsub("_+", "_", safe_name)
      safe_name <- gsub("^_|_\$", "", safe_name)
      safe_name
    }

    count_table <- read_delimited(count_matrix_file)

    if (ncol(count_table) < 3) {
      stop("Count matrix must contain a gene ID column and at least two sample columns.")
    }

    gene_id_column <- names(count_table)[1]
    rownames(count_table) <- make.unique(as.character(count_table[[gene_id_column]]))
    count_data <- count_table[, -1, drop = FALSE]
    count_data[] <- lapply(count_data, function(values) as.integer(round(as.numeric(values))))

    if (anyNA(count_data)) {
      stop("Count matrix contains missing or non-numeric values.")
    }

    sample_metadata <- read_delimited(metadata_file)

    if (!"sample" %in% names(sample_metadata)) {
      stop("Metadata must include a 'sample' column.")
    }

    sample_metadata[["sample"]] <- as.character(sample_metadata[["sample"]])

    if (anyDuplicated(sample_metadata[["sample"]])) {
      stop("Metadata contains duplicate sample names.")
    }

    missing_samples <- setdiff(sample_metadata[["sample"]], colnames(count_data))

    if (length(missing_samples) > 0) {
      stop("Metadata samples missing from count matrix: ", paste(missing_samples, collapse = ", "))
    }

    count_data <- count_data[, sample_metadata[["sample"]], drop = FALSE]
    rownames(sample_metadata) <- sample_metadata[["sample"]]

    for (column_name in names(sample_metadata)) {
      if (is.character(sample_metadata[[column_name]])) {
        sample_metadata[[column_name]] <- factor(sample_metadata[[column_name]])
      }
    }

    design <- as.formula(design_formula)
    missing_design_variables <- setdiff(all.vars(design), names(sample_metadata))

    if (length(missing_design_variables) > 0) {
      stop("Design variables missing from metadata: ", paste(missing_design_variables, collapse = ", "))
    }

    for (contrast_parts in contrasts) {
      contrast_variable <- contrast_parts[[1]]
      numerator_level <- contrast_parts[[2]]
      denominator_level <- contrast_parts[[3]]

      if (!contrast_variable %in% names(sample_metadata)) {
        stop("Contrast variable is absent from metadata: ", contrast_variable)
      }

      if (!contrast_variable %in% all.vars(design)) {
        stop("Contrast variable must be present in the DESeq2 design formula: ", contrast_variable)
      }

      if (!is.factor(sample_metadata[[contrast_variable]])) {
        sample_metadata[[contrast_variable]] <- factor(sample_metadata[[contrast_variable]])
      }

      if (!all(c(numerator_level, denominator_level) %in% levels(sample_metadata[[contrast_variable]]))) {
        stop("Contrast levels are absent from metadata column ", contrast_variable)
      }

      sample_metadata[[contrast_variable]] <- relevel(
        sample_metadata[[contrast_variable]],
        ref = denominator_level
      )
    }

    dds <- DESeqDataSetFromMatrix(
      countData = as.matrix(count_data),
      colData = sample_metadata,
      design = design
    )
    dds <- dds[rowSums(counts(dds)) > 0, ]

    if (nrow(dds) == 0) {
      stop("No genes remain after filtering zero-count rows.")
    }

    dds <- DESeq(dds)

    for (contrast_parts in contrasts) {
      contrast_name <- make_safe_name(contrast_parts)
      result <- results(dds, contrast = contrast_parts)
      result_table <- as.data.frame(result)
      result_table <- data.frame(
        gene_id = rownames(result_table),
        result_table,
        check.names = FALSE
      )
      result_table <- result_table[order(result_table[["padj"]], na.last = TRUE), , drop = FALSE]

      write.table(
        result_table,
        file = paste0("deseq2_results_", contrast_name, ".tsv"),
        sep = "\\t",
        quote = FALSE,
        row.names = FALSE
      )

      sig <- !is.na(result_table\$padj) & result_table\$padj < 0.05
      sig_fc <- sig & abs(result_table\$log2FoldChange) >= 1

      summary_table <- data.frame(
        contrast = contrast_name,
        genes_tested = sum(!is.na(result_table\$pvalue)),
        significant_padj_0.05 = sum(sig),
        significant_padj_0.05_log2FC_1 = sum(sig_fc),
        upregulated = sum(sig_fc & result_table\$log2FoldChange > 0),
        downregulated = sum(sig_fc & result_table\$log2FoldChange < 0)
      )

      write.table(
        summary_table,
        file = paste0("deseq2_summary_", contrast_name, ".tsv"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )

      pdf(paste0("deseq2_ma_plot_", contrast_name, ".pdf"))
      plotMA(result, ylim = c(-5, 5), main = paste(contrast_parts, collapse = " "))
      dev.off()

      volcano_df <- result_table[
        !is.na(result_table[["log2FoldChange"]]) &
        !is.na(result_table[["padj"]]),
        ,
        drop = FALSE
      ]

      volcano_df[["neglog10_padj"]] <- -log10(
        pmax(volcano_df[["padj"]], .Machine\$double.xmin)
      )

      volcano_df[["significant"]] <-
        volcano_df[["padj"]] < 0.05 &
        abs(volcano_df[["log2FoldChange"]]) >= 1

      # Squish extreme log2FC values to -10 / +10 for plotting
      volcano_df[["plot_log2FC"]] <- pmax(
        pmin(volcano_df[["log2FoldChange"]], 10),
        -10
      )

      pdf(paste0("deseq2_volcano_", contrast_name, ".pdf"))

      plot(
        volcano_df[["plot_log2FC"]],
        volcano_df[["neglog10_padj"]],
        pch = 16,
        cex = 0.5,
        col = ifelse(volcano_df[["significant"]], "red", "grey60"),
        xlab = "log2 Fold Change",
        ylab = "-log10 adjusted p-value",
        main = paste(contrast_parts, collapse = " "),
        xlim = c(-10, 10)
      )

      abline(h = -log10(0.05), lty = 2)
      abline(v = c(-1, 1), lty = 2)

      dev.off()

      png(
        paste0("deseq2_volcano_", contrast_name, ".png"),
        width = 1800,
        height = 1600,
        res = 200
      )

      plot(
        volcano_df[["plot_log2FC"]],
        volcano_df[["neglog10_padj"]],
        pch = 16,
        cex = 0.5,
        col = ifelse(volcano_df[["significant"]], "red", "grey60"),
        xlab = "log2 Fold Change",
        ylab = "-log10 adjusted p-value",
        main = paste(contrast_parts, collapse = " "),
        xlim = c(-10, 10)
      )

      abline(h = -log10(0.05), lty = 2)
      abline(v = c(-1, 1), lty = 2)

      dev.off()

    }

    normalized_counts <- counts(dds, normalized = TRUE)
    normalized_counts <- data.frame(
      gene_id = rownames(normalized_counts),
      normalized_counts,
      check.names = FALSE
    )

    write.table(
      normalized_counts,
      file = "deseq2_normalized_counts.tsv",
      sep = "\\t",
      quote = FALSE,
      row.names = FALSE
    )

    vst <- varianceStabilizingTransformation(dds, blind = FALSE)
    vst_counts <- assay(vst)
    vst_counts <- data.frame(
      gene_id = rownames(vst_counts),
      vst_counts,
      check.names = FALSE
    )

    write.table(
      vst_counts,
      file = "deseq2_vst_counts.tsv",
      sep = "\\t",
      quote = FALSE,
      row.names = FALSE
    )

    pca_groups <- unique(vapply(contrasts, function(contrast_parts) contrast_parts[[1]], character(1)))
    pdf("deseq2_pca.pdf")
    print(plotPCA(vst, intgroup = pca_groups))
    dev.off()

    # PCA with treatment ellipses + centroid distances
    group_var <- pca_groups[[1]]

    # Use the exact same PCA coordinates as DESeq2 plotPCA()
    pca_data <- plotPCA(
      vst,
      intgroup = group_var,
      returnData = TRUE
    )

    percent_var <- 100 * attr(pca_data, "percentVar")

    pca_df <- data.frame(
      sample = rownames(pca_data),
      PC1 = pca_data[["PC1"]],
      PC2 = pca_data[["PC2"]],
      group = pca_data[[group_var]],
      check.names = FALSE
    )

    pca_df[["group"]] <- factor(pca_df[["group"]])

    pca_df[["group"]] <- factor(pca_df[["group"]])


    # Group centroids
    centroids <- aggregate(
      cbind(PC1, PC2) ~ group,
      data = pca_df, FUN = mean
    )

    write.table(centroids,
      file = "deseq2_pca_centroids.tsv",
      sep = "\t", quote = FALSE, row.names = FALSE
    )


    # Euclidean distance between group centroids
    centroid_matrix <- as.matrix(
      centroids[, c("PC1", "PC2"), drop = FALSE]
    )

    rownames(centroid_matrix) <- as.character(centroids[["group"]])

    centroid_distances <- as.matrix(
      dist(centroid_matrix, method = "euclidean")
    )

    distance_table <- as.data.frame(
      as.table(centroid_distances),
      stringsAsFactors = FALSE
    )

    names(distance_table) <- c(
      "group1",
      "group2",
      "euclidean_distance_PC1_PC2"
    )

    # Keep each pair only once and remove self comparisons
    group_order <- as.character(centroids[["group"]])

    distance_table[["group1_index"]] <- match(
      distance_table[["group1"]],
      group_order
    )

    distance_table[["group2_index"]] <- match(
      distance_table[["group2"]],
      group_order
    )

    distance_table <- distance_table[
      distance_table[["group1_index"]] <
      distance_table[["group2_index"]],
      ,
      drop = FALSE
    ]

    distance_table <- distance_table[
      ,
      c(
        "group1",
        "group2",
        "euclidean_distance_PC1_PC2"
      ),
      drop = FALSE
    ]

    distance_table <- distance_table[
      order(
        distance_table[["euclidean_distance_PC1_PC2"]],
        decreasing = TRUE
      ),
      ,
      drop = FALSE
    ]

    write.table(distance_table,
      file = "deseq2_pca_centroid_distances.tsv",
      sep = "\t", quote = FALSE,  row.names = FALSE
    )


    # Choose ellipse vs convex hull based on group size
    group_sizes <- table(pca_df[["group"]])

    ellipse_groups <- names(group_sizes[group_sizes >= 4])
    hull_groups <- names(group_sizes[group_sizes >= 3 & group_sizes < 4])

    pca_ellipse_df <- pca_df[
      pca_df[["group"]] %in% ellipse_groups,
      ,
      drop = FALSE
    ]

    pca_hull_df <- pca_df[
      pca_df[["group"]] %in% hull_groups,
      ,
      drop = FALSE
    ]

    if (nrow(pca_hull_df) > 0) {
      pca_hull_df <- do.call(
        rbind,
        lapply(
          split(pca_hull_df, pca_hull_df[["group"]]),
          function(df) {
            df[chull(df[["PC1"]], df[["PC2"]]), , drop = FALSE]
          }
        )
      )
    }


    p_ellipse <- ggplot(
      pca_df,
      aes(
        x = PC1,
        y = PC2,
        color = group
      )
    ) +

      # Convex hull for small groups
      {
        if (nrow(pca_hull_df) > 0)
          geom_polygon(
            data = pca_hull_df,
            aes(
              x = PC1,
              y = PC2,
              group = group,
              fill = group
            ),
            alpha = 0.12,
            color = NA,
            inherit.aes = FALSE
          )
      } +

      # Ellipse for groups with enough samples
      {
        if (nrow(pca_ellipse_df) > 0)
          stat_ellipse(
            data = pca_ellipse_df,
            aes(
              x = PC1,
              y = PC2,
              group = group,
              color = group
            ),
            type = "norm",
            linewidth = 0.8,
            inherit.aes = FALSE
          )
      } +

      # Individual samples
      geom_point(size = 3) +

      # Group centroids
      geom_point(
        data = centroids,
        aes(
          x = PC1,
          y = PC2,
          color = group
        ),
        shape = 4,
        size = 5,
        stroke = 1.5,
        inherit.aes = FALSE
      ) +

      labs(
        x = paste0(
          "PC1: ",
          round(percent_var[[1]], 1),
          "% variance"
        ),
        y = paste0(
          "PC2: ",
          round(percent_var[[2]], 1),
          "% variance"
        ),
        color = group_var,
        fill = group_var,
        title = "PCA with group spread and centroids"
      ) +

      theme_bw()

    ggsave("deseq2_pca_ellipse.pdf",
      p_ellipse, width = 8, height = 6)

    capture.output(sessionInfo(), file = "deseq2_session_info.txt")
    RSCRIPT
    """
}



workflow deseq2_from_featurecounts {
    take:
    featurecounts_files
    metadata
    design_formula
    contrast

    main:
    merge_featurecounts(featurecounts_files)
    deseq2(merge_featurecounts.out.count_matrix, metadata, design_formula, contrast)

    emit:
    count_matrix = merge_featurecounts.out.count_matrix
    summary = deseq2.out.summary
    results = deseq2.out.results
    normalized_counts = deseq2.out.normalized_counts
    vst_counts = deseq2.out.vst_counts
    pca = deseq2.out.pca
    pca_ellipse = deseq2.out.pca_ellipse
    pca_distances = deseq2.out.pca_distances
    pca_centroids = deseq2.out.pca_centroids
    ma_plot = deseq2.out.ma_plot
    volcano = deseq2.out.volcano
    volcano_png = deseq2.out.volcano_png
    session_info = deseq2.out.session_info
}
