library(tidyverse)
## Correct for 3' bias by modelling coverage depth as a function
## of the distance from the 3' end, for each sample

file <- snakemake@input$cov
sample_id = snakemake@wildcards$sample_id

results <- read_tsv(file) |>
mutate(log_cov = log(cov+1), sample_id = sample_id) |>
group_by(gene) |>
mutate(dist = n() - pos + 1) |>
ungroup()

train <- results |>
group_by(gene) |>
# Trim 5' and 3' edges since coverage is often funky there
filter(pos >= min(pos) + min(100, n()/5) , pos <= max(pos) - min(100, n()/5)) |>
ungroup()

# Model expression across all genes as a function of
# distance from the 3' end, on a log scale (so exponential fit)
model <- lm(data = train, log_cov ~ dist)

results$predicted <- predict(model, newdata = results)

results <- results |>
mutate(predicted = exp(predicted)) |>
group_by(gene) |>
# Correct the data based off fit values
mutate(corrected = cov * mean(cov)/predicted)

results |> write_tsv(snakemake@output$corrected)
