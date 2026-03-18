library(readr)
library(dplyr)
library(magrittr)
library(stringr)
library(tibble)

# The below function is coded by Claude Opus 4.6
# SEE ALSO: https://link.springer.com/article/10.1007/BF01253782
#' Hierarchical Shapley for factorial ANOVA.
#'
#' Computes Shapley-style importance for factors with a permission structure
#' (some factors must "arrive" before others) and an optional stratification
#' variable to avoid fitting huge models.
#'
#' @param data       Data frame.
#' @param resp       Name of response variable (string).
#' @param players    Character vector of factor names (the Shapley players).
#' @param prereqs    Named list specifying permission structure. Each entry
#'                   maps a player to the players that must precede it.
#'                   E.g. list(sample = c("batch")) means sample can only
#'                   arrive after batch. Players not listed have no prereqs.
#' @param strat      Optional: name of the stratification variable (string).
#'                   Must also be one of the players. Models containing this
#'                   player are fit within each level and SS summed, avoiding
#'                   a giant design matrix.
#' @param model_for_subset  Function(player_names, data, resp) -> SS_reg.
#'                   If NULL, default builds full factorial among the factors.
#'                   Override for custom model structures.
#' @return List with shapley values, r2_full, and diagnostics.
hierarchical_shapley <- function(data, resp, players, prereqs = list(),
                                 strat = NULL) {

  p <- length(players)
  stopifnot(p <= 8)  # sanity: 8! = 40320 permutations

  # --- Permission structure: enumerate admissible orderings ---
  # An ordering is admissible if for every player, all its prereqs appear earlier.
  all_perms <- function(n) {
    if (n == 1) return(matrix(1, 1, 1))
    prev <- all_perms(n - 1)
    do.call(rbind, lapply(1:n, function(pos) {
      cbind(
        prev[,seq_len(pos-1),drop=FALSE],
        n,
        prev[,seq_len(ncol(prev)-pos+1) + pos -1, drop=FALSE]
      )
    }))
  }

  perms <- all_perms(p)

  # Filter to admissible orderings
  is_admissible <- function(perm) {
    pos_of <- setNames(seq_along(perm), players[perm])
    for (pl in names(prereqs)) {
      for (pr in prereqs[[pl]]) {
        if (pos_of[pl] <= pos_of[pr]) return(FALSE)
      }
    }
    TRUE
  }
  admissible <- which(apply(perms, 1, is_admissible))
  perms <- perms[admissible, , drop = FALSE]
  n_perms <- nrow(perms)
  cat(sprintf("  %d / %d orderings admissible\n", n_perms, factorial(p)))

  # --- Precompute R²(S) for all needed subsets ---
  # With permission structure, not all subsets appear in admissible orderings.
  # But with <=8 players, just precompute all 2^p - 1 subsets.

  grand_mean <- mean(data[[resp]])
  ss_total <- sum((data[[resp]] - grand_mean)^2)

  # Stratification setup
  if (!is.null(strat)) {
    stopifnot(strat %in% players)
    strat_data <- split(data, data[[strat]])
    ss_between <- ss_total - sum(sapply(strat_data, function(d)
      sum((d[[resp]] - mean(d[[resp]]))^2)))
  }

  # For a subset S of players, determine the effective factors in the model.
  # If strat is in S, all other players in S are used as factorial terms
  # fit within each stratum. If strat is not in S, fit on pooled data.
  compute_r2 <- function(pnames) {
    if (length(pnames) == 0) return(0)

    has_strat <- !is.null(strat) && (strat %in% pnames)
    model_factors <- setdiff(pnames, strat)

    if (has_strat && length(model_factors) == 0) {
      # Strat only: R² = SS_between / SS_total
      return(ss_between / ss_total)
    }

    f <- reformulate(paste(model_factors, collapse = " * "), resp)

    if (has_strat) {
      ss_reg_within <- sum(sapply(strat_data, function(d) {
        fit <- lm(f, data = d)
        sum((fitted(fit) - mean(d[[resp]]))^2)
      }))
      return((ss_between + ss_reg_within) / ss_total)
    }

    # No strat: pooled fit
    fit <- lm(f, data = data)
    sum((fitted(fit) - grand_mean)^2) / ss_total
  }

  # Make a lookup cache
  r2_cache <- new.env(hash = TRUE, parent = emptyenv())
  lookup <- function(idx) {
    if (length(idx) == 0) return(0)
    key <- paste(sort(idx), collapse = ",")
    if (is.null(r2_cache[[key]])) {
        r2_cache[[key]] <- compute_r2(players[idx])
    }
    r2_cache[[key]]
  }

  # --- Shapley values via averaging over admissible orderings ---
  shapley <- setNames(numeric(p), players)
  for (r in 1:n_perms) {
    perm <- perms[r, ]
    for (pos in 1:p) {
      j <- perm[pos]
      before <- if (pos == 1) integer(0) else perm[1:(pos-1)]
      marginal <- lookup(c(before, j)) - lookup(before)
      shapley[players[j]] <- shapley[players[j]] + marginal / n_perms
    }
  }

  list(
    shapley = shapley,
    r2_full = unname(lookup(1:p)),
    n_admissible = n_perms,
    r2_subsets = r2_cache %>% as.list()
  )
}


cov <- read_tsv("results/transcript_coverage/liver.transcript_coverage.txt.gz")
sample_info <- read_tsv("results/liver.sample_info.txt")

cov_downsampled <- cov %>%
    group_by(gene, sample) %>%
    mutate( mean_gene_cov = mean(cov) ) %>%
    filter( pos %% 100 == 50 ) %>%
    mutate( pos_gene = paste0(gene, pos)) %>%
    mutate( cov_normalized = cov / mean_gene_cov) %>%
    left_join(
        sample_info %>% dplyr::select(sample = ID, study),
        by = "sample",
    )

# Shapley values
res <- hierarchical_shapley(
    cov_downsampled,
    "cov_normalized",
    c("study", "sample", "pos_gene"),
    strat = "pos_gene",
    prereqs = list(
        study = c("pos_gene"),
        sample = c("study", "pos_gene"),
        pos_gene = c()
    )
)
tibble(variable = names(res$shapley), shapley_value=res$shapley) %>% write_tsv("results/scratch/coverage.GEO.shapley_values.txt")


## SAME BUT FOR THE SEQC / UHR DATASET
cov <- read_tsv("results/transcript_coverage/UHR.transcript_coverage.txt.gz")
sample_info <- read_tsv("results/UHR.sample_info.txt") %>%
    mutate(
        sequencing_batch = case_when(
            library_id %in% c(1,2,3,4) ~ str_split_i(study, "_", 2),
            ID == "SRX302574" ~ "CNL",
            ID == "SRX302194" ~ "BGI",
            ID == "SRX302938" ~ "MAY",
        ),
        library_prep_batch = str_split_i(study, "_", 2),
    )

cov_downsampled <- cov %>%
    group_by(gene, sample) %>%
    mutate( mean_gene_cov = mean(cov) ) %>%
    filter( pos %% 100 == 50 ) %>%
    mutate( cov_normalized = cov / mean_gene_cov) %>%
    mutate( pos_gene = paste0(gene, pos)) %>%
    left_join(
        sample_info %>% dplyr::select(sample = ID, study, library_prep_batch, sequencing_batch),
        by = "sample",
    )

res <- hierarchical_shapley(
    cov_downsampled,
    "cov_normalized",
    c("library_prep_batch", "sequencing_batch", "sample", "pos_gene"),
    strat = "pos_gene",
    prereqs = list(
        sample = c("library_prep_batch", "sequencing_batch"),
        library_prep_batch = c("pos_gene"),
        sequencing_batch = c("pos_gene"),
        pos_gene = c()
    )
)
tibble(variable = names(res$shapley), shapley_value=res$shapley) %>% write_tsv("results/scratch/coverage.UHR.shapley_values.txt")

## SEQC A + B, UHR + HBRR comparison
cov1 <- read_tsv("results/transcript_coverage/UHR.transcript_coverage.txt.gz")
cov2 <- read_tsv("results/transcript_coverage/HBRR.transcript_coverage.txt.gz")
in_both <- intersect(cov1$gene %>% unique(), cov2$gene %>% unique())
cov <- bind_rows(
    cov1 %>% filter(gene %in% in_both),
    cov2 %>% filter(gene %in% in_both),
)
sample_info <- bind_rows(
        read_tsv("results/UHR.sample_info.txt"),
        read_tsv("results/HBRR.sample_info.txt"),
    )%>%
    mutate(
        sequencing_batch = case_when(
            ID == "SRX302574" ~ "CNL",
            ID == "SRX302194" ~ "BGI",
            ID == "SRX302938" ~ "MAY",
            ID == 'SRX302288' ~ "BGI",
            ID == 'SRX303021' ~ "MAY",
            ID == 'SRX302651' ~ "CNL",
            TRUE ~ str_split_i(study, "_", 2),
        ),
        library_prep_batch = str_split_i(study, "_", 2),
    )

cov_downsampled <- cov %>%
    group_by(gene, sample) %>%
    mutate( mean_gene_cov = mean(cov) ) %>%
    filter( pos %% 100 == 50 ) %>%
    mutate( cov_normalized = cov / mean_gene_cov) %>%
    mutate( pos_gene = paste0(gene, pos)) %>%
    left_join(
        sample_info %>% select(sample = ID, study, library_prep_batch, sequencing_batch, tissue),
        by = "sample",
    )

res <- hierarchical_shapley(
    cov_downsampled,
    "cov_normalized",
    c("library_prep_batch", "sequencing_batch", "sample", "tissue", "pos_gene"),
    strat = "pos_gene",
    prereqs = list(
        sample = c("library_prep_batch", "sequencing_batch", 'tissue'),
        library_prep_batch = c("pos_gene"),
        sequencing_batch = c("pos_gene"),
        tissue = c("pos_gene"),
        pos_gene = c()
    )
)
tibble(variable = names(res$shapley), shapley_value=res$shapley) %>% write_tsv("results/scratch/coverage.UHR_HBRR.shapley_values.txt")
