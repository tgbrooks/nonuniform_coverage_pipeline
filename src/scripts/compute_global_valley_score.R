library(tidyverse)

KERNEL_SIZE <- 101
CONV_SD <- 10

cov_file <- snakemake@input$cov
#cov_file <- "results/alpine_fits/liver/SRX11694499.coverage_table.txt"

### BEGIN Preprocessing data
coverage <- read_tsv(cov_file) |>
    rename(gene_id = gene, cov = actual)

if ("sample" %in% colnames(coverage)) {
    group_vars = c("sample", "gene_id")
} else {
    group_vars = c("gene_id")
}

coverage_stats <- coverage |>
    group_by(across(all_of(group_vars))) |>
    summarize(
        transcript_length = max(pos)+1,
        mean_cov = mean(cov),
    )

gene_ids <- unique(coverage$gene_id)

### END Preprocessing data

### BEGIN Functions
## Convolution functions
create_gaussian_kernel <- function(conv_sd, kernel_size) {
  x <- seq(-kernel_size / 2, kernel_size / 2, length.out = kernel_size)
  kernel <- exp(-0.5 * (x / conv_sd)^2)
  kernel <- kernel / sum(kernel)  # Normalize the kernel
  return(kernel)
}

get_conv_result <- function(x, conv_sd=sd) {
  conv_result <- convolve(x, create_gaussian_kernel(conv_sd, KERNEL_SIZE), type = "open")
  result_length <- length(conv_result)
  rem <- floor(KERNEL_SIZE/2)
  output <- conv_result[(rem + 1) : (result_length - rem)]
  return(output)
}

## Two-pointer Valley and Peak finders
tp_valley_finder <- function(data) {

  len_gene <- length(data)
  data <- c(0,0,data,0,0)

  result <- tibble(left = numeric(), right = numeric(), conv_depth = numeric())
  left = 1
  #min_valley_depth <- sqrt(mean(data))
  min_valley_depth <- mean(data)/10
  first_peak_found <- FALSE
  running_min <- data[left]

  for (right in 2:(length(data)-1)) {
    running_min <- min(running_min, data[right])

    ## Out of valley check
    if((data[right - 1] < data[right]) && (data[right] >= data[right + 1])){
      score <- min(data[left]-running_min, data[right]-running_min)
      #cat("Left:", left, ", Right:", right, ", Running Min:", running_min, ", Score:", score, "\n")
      if (score > min_valley_depth) {
        result <- add_row(result, left = left, right = right, conv_depth = score)
        left <- right
        running_min <- data[left]
      }
      else if( (score <= min_valley_depth) && (data[left] < data[right]) ){
        left <- right
        running_min <- data[left]
      }
    }

    if ((data[right-1] <= data[right]) && (data[right] > data[right+1])){
      if (!first_peak_found) {
        first_peak_found <- TRUE
        left <- right
        running_min <- data[left]
      }
    }
  }

  result <- result |>
    mutate(left = pmin(pmax(left - 2,1), len_gene), right = pmax(pmin(right - 2, len_gene),1))

  return(result)
}

tp_peak_finder <- function(data) {
  M <- max(data)
  result <- tp_valley_finder(M - data)
  return(result)
}

## Helper function for finding intervals around valleys and peaks: from a vector of Boolean, extracts all the intervals of consecutive "TRUE"s
get_interval_list <- function(mask){
  interval_list <- list()
  k <- 1
  l <- 1
  in_interval <- FALSE
  if(mask[l] == TRUE){
    in_interval <- TRUE
  }else{
    in_interval <- FALSE
  }

  mask <- c(mask, FALSE)

  for(r in 1:length(mask)){
    if(in_interval && (mask[r] == FALSE)){
      interval_list[[k]] <- c(l,r-1)
      k <- k+1
      in_interval <- FALSE
    }
    else if(!in_interval && (mask[r] == TRUE)){
      l <- r
      in_interval <- TRUE
    }
  }
  return(interval_list)
}


#### Run smoothing convolution

cov_conv <- coverage |>
    group_by(gene_id) |>
    mutate(conv = get_conv_result(cov, conv_sd = CONV_SD))

## Get all valleys (with Gauss smoothing)
valley_results <- list()
valley_results <- cov_conv |>
    group_by(across(all_of(group_vars))) |>
    group_modify(function(data, groups) { 
        tp_valley_finder(data$conv)
    })

all_valley_results <- valley_results |>
  select(all_of(group_vars), left, right) |>
  mutate(optimum_type = "Valley")

#### Find local regions around Peaks and Valleys

# NOTE: this is not used in this script, which is just for computing the global valley score

#### Get intervals for valleys
##temp <- list()
##for(i in 1:nrow(all_valley_results)){
##  g_id <- all_valley_results$gene_id[i]
##  left <- all_valley_results$left[i]
##  right <- all_valley_results$right[i]
##
##  gene_data <- (filter(cov_conv, gene_id == g_id)[["conv"]])
##  left_end <- gene_data[left]
##  right_end <- gene_data[right]
##
##  m <- min(gene_data[left:right])
##  low_thresh <- m + min(left_end - m, right_end - m)/3
##
##  ## Find all the sub-intervals satisfying x in SI <==> x <= thresh
##  low_mask <- (gene_data[left:right] < low_thresh)
##  list_of_intervals <- get_interval_list(low_mask)
##
##  for(j in 1:length(list_of_intervals)){
##    left_end_pt <- list_of_intervals[[j]][1] + left - 1
##    right_end_pt <- list_of_intervals[[j]][2] + left - 1
##    temp[[length(temp)+1]] <- tibble(gene_id = g_id, left = left_end_pt, right = right_end_pt, type = "Valley")
##  }
##}
##valley_sub_intervals <- bind_rows(temp)
##
#### Get intervals for peaks
##temp <- list()
##for(i in 1:nrow(all_peak_results)){
##  g_id <- all_peak_results$gene_id[i]
##  left <- all_peak_results$left[i]
##  right <- all_peak_results$right[i]
##
##  gene_data <- (filter(cov_conv, gene_id == g_id)[["conv"]])
##  left_end <- gene_data[left]
##  right_end <- gene_data[right]
##
##  M <- max(gene_data[left:right])
##  high_thresh <- M - max(M - left_end, M - right_end)/3
##
##  ## Find all the sub-intervals satisfying x in SI <==> x >= thresh
##  high_mask <- (gene_data[left:right] > high_thresh)
##  list_of_intervals <- get_interval_list(high_mask)
##
##  #temp <- tibble()
##  for(j in 1:length(list_of_intervals)){
##    left_end_pt <- list_of_intervals[[j]][1] + left - 1
##    right_end_pt <- list_of_intervals[[j]][2] + left - 1
##    temp[[length(temp)+1]] <- tibble(gene_id = g_id, left = left_end_pt, right = right_end_pt, type = "Peak")
##  }
##}
##peaks_sub_intervals <- bind_rows(temp)
##
#### Combine the valley and peak intervals data
##
##all_peaks_valleys_intervals <- bind_rows(peaks_sub_intervals, valley_sub_intervals) |>
##  arrange(gene_id)
##
###write_csv(all_peaks_valleys_intervals, "all_peaks_valleys_intervals.csv")

### COMPUTE VALLEY SCORES

valley_integral <- function(left, right, data){
  m <- min(data[left], data[right])
  temp <- data[left:right]
  return(sum(pmax(m - temp, 0)))
}


#local_valley_score <- (valley_integral(left, right, conv_data) / mean(data)) / transcript_length
cov_by_gene <- list()
for (g_id in gene_ids) {
    cov_by_gene[[g_id]] <- (coverage |> filter(gene_id == g_id))[['cov']]
}

local_valley_score <- all_valley_results |>
    rowwise() |>
    mutate(
        lvs = valley_integral(
            left,
            right,
            cov_by_gene[[gene_id]]
        )
    )

global_valley_score <- local_valley_score |>
    left_join(coverage_stats, by=group_vars) |>
    group_by(across(all_of(group_vars))) |>
    summarize(
        global_valley_score = sum(lvs / mean_cov / transcript_length),
    )

if ("sample" %in% colnames(coverage)) {
    results <- global_valley_score |>
        group_by(sample) |>
        summarize(global_valley_score = sum(global_valley_score))
} else {
    results <- tibble(
        global_valley_score = sum(global_valley_score$global_valley_score),
    )
}
write_tsv(results, snakemake@output$gvs)
