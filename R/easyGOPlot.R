#' Plot GO enrichment results
#'
#' Plots the most significant GO terms from [easytopGO()] or a VEuPathDB
#' enrichment table. Input p-values are sorted before selecting the top terms.
#' Redundant rows are identified using GO.ID when available, otherwise the term
#' description; this is a display-only reduction, not a change to enrichment.
#'
#' @param res A data.frame of GO enrichment results.
#' @param title Plot title.
#' @param limit Maximum number of terms to display (positive integer).
#' @param desc Name of the term-description column (e.g. "Term" or "Name").
#' @param genecounts Name of the gene-count column (e.g. "Significant" or
#'   "Result count"). For KS, Significant may refer to a nominal gene cutoff.
#' @param sortby Name of a numeric GO p-value column. When NULL (default), uses
#'   the first available among padj, fisher, ks, P-value and p-value.
#' @param min_genes Optional minimum count for plotting (NULL means no filter).
#'   Applied only to the displayed terms, not before GO enrichment/FDR testing.
#'
#' @return A ggplot2 object. Empty results produce a labelled empty plot.
#' @importFrom rlang .data
#' @export
#' @examples
#' \dontrun{
#' gores <- easytopGO(
#'   geneID = down_padj, bkggset = background,
#'   gaf = "annotation.gaf.gz", useBiomart = FALSE,
#'   useGAF = TRUE, stats = "fisher"
#' )
#' easyGOPlot(gores, sortby = "fisher", min_genes = 3)
#' # For KS output:
#' # easyGOPlot(ks_result, sortby = "ks")
#' # For a VEuPathDB result with numeric p-values:
#' # easyGOPlot(veupath, desc = "Name", genecounts = "Result count",
#' #            sortby = "P-value")
#' }
easyGOPlot <- function(res, title = "GO Biological processes", limit = 20L,
                       desc = "Term", genecounts = "Significant",
                       sortby = NULL, min_genes = NULL) {
  if (!is.data.frame(res)) stop("res must be a data.frame or tibble")
  if (anyDuplicated(names(res))) stop("res has duplicate column names")
  if (!is.numeric(limit) || length(limit) != 1L || is.na(limit) ||
    !is.finite(limit) || limit < 1 || limit != floor(limit)) {
    stop("limit must be a positive whole number")
  }
  if (!is.null(min_genes) &&
    (!is.numeric(min_genes) || length(min_genes) != 1L ||
      is.na(min_genes) || !is.finite(min_genes) || min_genes < 0 ||
      min_genes != floor(min_genes))) {
    stop("min_genes must be NULL or a non-negative whole number")
  }
  if (!is.character(desc) || length(desc) != 1L || is.na(desc) ||
    !nzchar(desc) || !is.character(genecounts) ||
    length(genecounts) != 1L || is.na(genecounts) || !nzchar(genecounts)) {
    stop("desc and genecounts must each be one column name")
  }

  # Pick the most appropriate supplied GO p-value column automatically.
  if (is.null(sortby)) {
    options <- c("padj", "fisher", "ks", "P-value", "p-value", "P value")
    found <- intersect(options, names(res))
    if (!length(found)) {
      stop("Cannot detect a GO p-value column; supply sortby explicitly")
    }
    sortby <- found[1L]
  }
  if (!is.character(sortby) || length(sortby) != 1L || is.na(sortby) ||
    !nzchar(sortby)) {
    stop("sortby must be one p-value column name")
  }
  needed <- c(desc, genecounts, sortby)
  if (!all(needed %in% names(res))) {
    stop("Missing column(s): ", paste(setdiff(needed, names(res)), collapse = ", "))
  }

  empty_plot <- function(note) {
    ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = note, size = 4) +
      ggplot2::theme_void(base_size = 12) +
      ggplot2::labs(title = title)
  }
  if (!nrow(res)) {
    return(empty_plot("No enriched GO terms to display"))
  }
  res <- as.data.frame(res)

  # Work on copies: avoid overwriting the raw enrichment p-value column.
  p_text <- trimws(as.character(res[[sortby]]))
  less_than <- !is.na(p_text) & grepl("^<\\s*", p_text)
  p_text[less_than] <- sub("^<\\s*", "", p_text[less_than])
  pval <- suppressWarnings(as.numeric(p_text))
  malformed <- !is.na(p_text) & nzchar(p_text) & is.na(pval)
  if (any(malformed)) {
    stop(
      "Non-numeric GO p-value(s) in '", sortby, "': ",
      paste(utils::head(unique(p_text[malformed]), 3L), collapse = ", ")
    )
  }
  if (any(less_than)) {
    warning("GO p-values printed as '< x' are plotted at their stated ",
      "upper bound x; exact values are unavailable.",
      call. = FALSE
    )
  }
  if (any(!is.na(pval) & (!is.finite(pval) | pval < 0 | pval > 1))) {
    stop("GO p-values must be finite numbers between 0 and 1")
  }

  count_text <- gsub(",", "", trimws(as.character(res[[genecounts]])),
    fixed = TRUE
  )
  counts <- suppressWarnings(as.numeric(count_text))
  malformed_count <- !is.na(count_text) & nzchar(count_text) & is.na(counts)
  if (any(malformed_count) ||
    any(!is.na(counts) & (!is.finite(counts) | counts < 0))) {
    stop("'", genecounts, "' must contain finite, non-negative numeric counts")
  }
  term <- trimws(as.character(res[[desc]]))
  keep <- !is.na(pval) & !is.na(counts) & !is.na(term) & nzchar(term)
  if (!is.null(min_genes)) keep <- keep & counts >= min_genes
  if (any(!keep)) {
    message(sum(!keep), " term(s) excluded from plot (missing values or count filter)")
  }
  res <- res[keep, , drop = FALSE]
  pval <- pval[keep]
  counts <- counts[keep]
  term <- term[keep]
  if (!nrow(res)) {
    return(empty_plot("No GO terms pass the plotting filters"))
  }

  # Crucial: rank terms BEFORE choosing top 'limit' entries.
  ord <- order(pval, -counts, term)
  res <- res[ord, , drop = FALSE]
  pval <- pval[ord]
  counts <- counts[ord]
  term <- term[ord]

  # A GO identifier is a safer duplicate key than description alone.
  # VEuPathDB tables without GO.ID are deduplicated by description.
  if ("GO.ID" %in% names(res)) {
    go_id <- trimws(as.character(res[["GO.ID"]]))
    key <- ifelse(!is.na(go_id) & nzchar(go_id),
      paste0("ID:", go_id), paste0("TERM:", term)
    )
  } else {
    key <- paste0("TERM:", term)
  }
  unique_row <- !duplicated(key)
  res <- res[unique_row, , drop = FALSE]
  pval <- pval[unique_row]
  counts <- counts[unique_row]
  term <- term[unique_row]

  nplot <- min(nrow(res), as.integer(limit))
  res <- res[seq_len(nplot), , drop = FALSE]
  pval <- pval[seq_len(nplot)]
  counts <- counts[seq_len(nplot)]
  term <- term[seq_len(nplot)]

  # Make identical descriptions for distinct GO IDs unambiguous on the plot.
  label <- term
  if (anyDuplicated(label) && "GO.ID" %in% names(res)) {
    duplicated_labels <- duplicated(label) | duplicated(label, fromLast = TRUE)
    label[duplicated_labels] <- paste0(
      label[duplicated_labels], " [", res[["GO.ID"]][duplicated_labels], "]"
    )
  }
  label <- make.unique(label)

  # Zero can arise from numeric underflow, and -log10(0) is infinite.
  # Cap the visual score only; do not modify the original results/p-values.
  zero_p <- any(pval == 0)
  res$.easygo_score <- -log10(pmax(pval, .Machine$double.xmin))
  res$.easygo_count <- counts
  # First selected row is most significant; place it at the top after flipping.
  res$.easygo_term <- factor(label, levels = rev(label))

  ranges <- unique(round(seq(min(counts), max(counts),
    length.out = min(3L, nplot)
  )))
  caption <- paste0("Reference lines: p = 0.05, 0.01, 0.001 (", sortby, ")")
  if (zero_p) {
    caption <- paste0(
      caption, "; p = 0 visually capped at .Machine$double.xmin"
    )
  }

  ggplot2::ggplot(res, ggplot2::aes(
    x = .data[[".easygo_term"]], y = .data[[".easygo_score"]],
    size = .data[[".easygo_count"]], fill = .data[[".easygo_score"]]
  )) +
    ggplot2::geom_hline(
      yintercept = -log10(0.05), linetype = "dotted",
      colour = "black", linewidth = 0.35, alpha = 0.5
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(0.01), linetype = "longdash",
      colour = "black", linewidth = 0.45, alpha = 0.5
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(0.001), linetype = "solid",
      colour = "black", linewidth = 0.55, alpha = 0.5
    ) +
    ggplot2::geom_point(shape = 21, colour = "black", stroke = 0.25) +
    ggplot2::scale_size_continuous(
      range = c(3, 8), name = "No. of Genes", breaks = ranges
    ) +
    ggplot2::scale_fill_gradient(
      low = "#7794F3", high = "#F63C3D",
      name = paste0("-log10(", sortby, ")")
    ) +
    ggplot2::scale_x_discrete(labels = scales::label_wrap(40)) +
    ggplot2::labs(
      x = NULL, y = paste0("-log10(", sortby, ")"),
      title = title,
      subtitle = paste("Top", nplot, "terms ordered by", sortby),
      caption = caption
    ) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      legend.position = "right",
      plot.title = ggplot2::element_text(size = 12, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 10),
      plot.caption = ggplot2::element_text(size = 9),
      axis.text.x = ggplot2::element_text(size = 10),
      axis.text.y = ggplot2::element_text(size = 10, face = "bold"),
      axis.title = ggplot2::element_text(size = 11, face = "bold"),
      axis.line = ggplot2::element_line(colour = "black"),
      legend.key = ggplot2::element_blank(),
      legend.key.size = grid::unit(0.7, "cm"),
      legend.text = ggplot2::element_text(size = 10)
    )
}
