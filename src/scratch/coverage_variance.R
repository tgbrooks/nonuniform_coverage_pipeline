library(readr)
library(dplyr)
library(magrittr)
library(stringr)
library(tibble)

# The below function is coded by Claude Opus 4.6
#' Grouped Shapley (LMG) for factorial ANOVA with a stratification variable.
#'
#' Players are the factors listed on the RHS of `formula` plus `group`.
#' For a subset S of players, the model includes all interactions among
#' the factors in S. When `group` is in S, we use the stratification trick
#' (fit within each group level, sum SS_reg) to avoid a giant design matrix.
#'
#' @param data   Data frame
#' @param formula  Formula like y ~ U * V * W (just the non-group factors)
#' @param group  Name of the grouping variable (string)
#' @return List with shapley values, R² of full model, per-subset R²
grouped_shapley <- function(data, formula, group) {
  resp <- all.vars(formula)[1]
  factors <- all.vars(formula)[-1]

  # Players: factors + group. Index group as last player.
  players <- c(factors, group)
  p <- length(players)
  gi <- p  # group index

  grp_data <- split(data, data[[group]])

  # SS_total: within-group (this is what stratified models explain)
  # Plus between-group (what group main effect explains)
  grand_mean <- mean(data[[resp]])
  ss_total <- sum((data[[resp]] - grand_mean)^2)
  ss_within <- sum(sapply(grp_data, function(d) sum((d[[resp]] - mean(d[[resp]]))^2)))
  ss_between <- ss_total - ss_within

  # Build interaction formula from a set of factor names
  make_formula <- function(fnames) {
    if (length(fnames) == 0) return(NULL)
    # Full factorial among fnames
    reformulate(paste(fnames, collapse = " * "), resp)
  }

  # Compute SS_reg for a subset of players (given as indices into `players`)
  compute_ss_reg <- function(idx) {
    if (length(idx) == 0) return(0)

    has_group <- gi %in% idx
    facs <- players[setdiff(idx, gi)]  # non-group factors in this subset

    if (has_group && length(facs) == 0) {
      # Only group: SS_reg = SS_between
      return(ss_between)
    }

    if (has_group) {
      # Stratified: fit y ~ facs (full factorial) within each group
      f <- make_formula(facs)
      ss_reg_within <- sum(sapply(grp_data, function(d) {
        fit <- lm(f, data = d)
        sum((fitted(fit) - mean(d[[resp]]))^2)
      }))
      # Group also explains between-group variance
      return(ss_between + ss_reg_within)
    }

    # No group: fit on pooled data
    f <- make_formula(facs)
    fit <- lm(f, data = data)
    sum((fitted(fit) - grand_mean)^2)
  }

  # Enumerate all non-empty subsets, compute R²
  subsets <- unlist(lapply(1:p, function(k) combn(p, k, simplify = FALSE)), recursive = FALSE)
  ss_reg_vals <- setNames(
    sapply(subsets, compute_ss_reg),
    sapply(subsets, function(s) paste(sort(s), collapse = ","))
  )
  r2_vals <- ss_reg_vals / ss_total

  # Shapley values
  lookup <- function(s) if (length(s) == 0) 0 else r2_vals[paste(sort(s), collapse = ",")]

  shapley <- setNames(numeric(p), players)
  for (j in 1:p) {
    others <- setdiff(1:p, j)
    val <- 0
    for (k in 0:length(others)) {
      S_list <- if (k == 0) list(integer(0)) else combn(others, k, simplify = FALSE)
      w <- 1 / (choose(p - 1, k) * p)
      for (S in S_list) val <- val + w * (lookup(c(S, j)) - lookup(S))
    }
    shapley[j] <- val
  }

  list(
    shapley = shapley,
    r2_full = unname(lookup(1:p)),
    r2_subsets = r2_vals
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
res <- grouped_shapley(cov_downsampled, cov_normalized ~ study * sample, "pos_gene")
tibble(variable = names(res$shapley), shapley_value=res$shapley) %>% write_tsv("results/scratch/coverage.GEO.shapley_values.txt")


## SAME BUT FOR THE SEQC / UHR DATASET
cov <- read_tsv("results/transcript_coverage/UHR.transcript_coverage.txt.gz")
sample_info <- read_tsv("results/UHR.sample_info.txt") %>%
    mutate(
        library_prep_batch = case_when(
            library_id %in% c(1,2,3,4) ~ str_split_i(study, "_", 2),
            ID == "SRX302574" ~ "CNL",
            ID == "SRX302194" ~ "BGI",
            ID == "SRX302938" ~ "MAY",
        ),
        sequencing_batch = str_split_i(study, "_", 2),
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

res <- grouped_shapley(cov_downsampled, cov_normalized ~ library_prep_batch + sequencing_batch + sample, "pos_gene")
tibble(variable = names(res$shapley), shapley_value=res$shapley) %>% write_tsv("results/scratch/coverage.UHR")

## SEQC A + B, UHR + HBRR comparison
cov <- bind_rows(
    read_tsv("results/transcript_coverage/UHR.transcript_coverage.txt.gz"),
    read_tsv("results/transcript_coverage/HBRR.transcript_coverage.txt.gz"),
)
sample_info <- bind_rows(
        read_tsv("results/UHR.sample_info.txt"),
        read_tsv("results/HBRR.sample_info.txt"),
    )%>%
    mutate(
        library_prep_batch = case_when(
            ID == "SRX302574" ~ "CNL",
            ID == "SRX302194" ~ "BGI",
            ID == "SRX302938" ~ "MAY",
            ID == 'SRX302288' ~ "BGI",
            ID == 'SRX303021' ~ "MAY",
            ID == 'SRX302651' ~ "CNL",
            TRUE ~ str_split_i(study, "_", 2),
        ),
        sequencing_batch = str_split_i(study, "_", 2),
    )

cov_downsampled <- cov %>%
    group_by(gene, sample) %>%
    mutate( mean_gene_cov = mean(cov) ) %>%
    filter( pos %% 100 == 50 ) %>%
    mutate( cov_normalized = cov / mean_gene_cov) %>%
    left_join(
        sample_info %>% select(sample = ID, study, library_prep_batch, sequencing_batch),
        by = "sample",
    )

# Position-wise fit
# These are done base-to-base since we want to account for positional variation
# but the full model is slow to run if we include a separate term for each base
# So instead we run one lm() for each base but explicitly subtract out the grand mean
# and add a dummy 'position' varaible that (since intercept is not included)
# accounts for the difference from the grand mean, thereby measuring variance of
# position.
grand_mean <- cov_downsampled$cov_normalized %>% mean()
models = list(
    base = cov_normalized - grand_mean ~ 0,
    position = cov_normalized - grand_mean ~ 1,
    library_prep = cov_normalized - grand_mean ~ library_prep_batch,
    sequencing_batch = cov_normalized - grand_mean ~ library_prep_batch + sequencing_batch,
    tissue = cov_normalized - grand_mean ~ library_prep_batch + sequencing_batch + tissue,
    tech_rep = cov_normalized - grand_mean ~ sample
)
model_fit <- function(data, keys) {
    fits <- lapply(models, function(model) { lm(model, data) })
    names(fits) <- NULL
    aov <- do.call(anova, fits)
    rownames(aov) <- names(models)
    tibble(
        gene_id = keys$gene,
        pos = keys$pos,
        variable = rownames(aov),
        sum_sq = aov$`Sum of Sq`,
    )
}
squared_sums <- cov_downsampled %>% group_by(gene, pos) %>% group_map(model_fit) %>% bind_rows()

overall <- squared_sums %>%
    filter(variable != 'base') %>%
    group_by(variable) %>%
    summarize(sum_sq = sum(sum_sq)) %>%
    mutate(pct_variance = 100*sum_sq / sum(sum_sq)) %>%
    arrange(factor(variable, levels=names(models)))


write_tsv(overall, "results/scratch/coverage.UHR_HBRR.percent_variance_explained.txt")
