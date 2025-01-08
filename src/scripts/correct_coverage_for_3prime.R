library(tidyverse)
## Correct for 3' bias by modelling coverage depth as a function
## of the distance from the 3' end, for each sample

s <- list()
for(sample_id in snakemake@params$sample_ids){
  file <- paste0("results/alpine_fits/", snakemake@wildcards$tissue, "/", sample_id, ".coverage_table.txt")
  temp <- read_tsv(file) |>
    mutate(log_actual = log(actual+1), sample_id = sample_id) |>
    group_by(gene) |>
    mutate(dist = n() - pos + 1) |>
    ungroup()

  train <- temp |>
    group_by(gene) |>
    # Trim 5' and 3' edges since coverage is often funky there
    filter(pos >= min(pos) + min(100, n()/5) , pos <= max(pos) - min(100, n()/5)) |>
    ungroup()

  # Model expression across all genes as a function of
  # distance from the 3' end, on a log scale (so exponential fit)
  model <- lm(data = train, log_actual ~ dist)

  temp$predicted <- predict(model, newdata = temp)

  temp <- temp |>
    mutate(predicted = exp(predicted)) |>
    group_by(gene) |>
    # Correct the data based off fit values
    mutate(corrected = actual * mean(actual)/predicted)

  s[[length(s) + 1]] <- temp
}

results <- bind_rows(s)

results |> write_tsv(snakemake@output$corrected)
