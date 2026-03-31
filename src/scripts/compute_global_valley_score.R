library(tidyverse)

KERNEL_SIZE <- 101
CONV_SD <- 10

cov_file <- snakemake@input$cov
#cov_file <- "results/alpine_fits/liver/SRX11694499.coverage_table.txt"

### BEGIN Preprocessing data
coverage <- read_tsv(cov_file) |>
    rename(gene_id = gene)

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

alt_valley_finder <- function(data) {
    len_gene <- length(data)
    data <- c(0,data,0)
    min_valley_depth <- mean(data)/10

    cum_left_max <- cummax(data)
    cum_right_max <- cummax(rev(data)) |> rev()
    water_height <- pmin(cum_right_max, cum_left_max)
    water_depth <- water_height - data # depth above the terrain
    y <- rle(water_depth > 0)
    valley_rights <- cumsum(y$lengths)[y$values == TRUE]
    valley_lefts <- valley_rights - y$lengths[y$values == TRUE]
    print(water_depth)
    print(y)

    valley_depths <- rep(0, length(valley_rights))
    valley_areas <- rep(0, length(valley_rights))
    for (i in seq_along(valley_rights)) {
        left <- valley_lefts[i]
        right <- valley_rights[i]
        valley_depths[i] <- max(water_depth[left:right])
        valley_areas[i] <- sum(water_depth[left:right])
    }

    results <- tibble(
        left = valley_lefts - 1, # minus the padding
        right = valley_rights, # inclusive of the endpoint
        depth = valley_depths,
        area = valley_areas,
    )
    return(results |> filter(depth > min_valley_depth))
}

## Two-pointer Valley and Peak finders
tp_valley_finder_orig <- function(data) {
  len_gene <- length(data)
  data <- c(0,0,data,0,0)

  result <- tibble(left = numeric(), right = numeric(), conv_depth = numeric())
  left = 1
  min_valley_depth <- mean(data)/10
  first_peak_found <- FALSE
  running_min <- data[left]

  for (right in 2:(length(data)-1)) {
    running_min <- min(running_min, data[right])

    # Out of valley check
    if((data[right - 1] < data[right]) && (data[right] >= data[right + 1])){
      score <- min(data[left]-running_min, data[right]-running_min)
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

tp_valley_finder <- function(data) {

    len_gene <- length(data)
    data <- c(0,0,data,0,0)

    result <- tibble(left = numeric(), right = numeric(), conv_depth = numeric())
    left = 1
    #min_valley_depth <- sqrt(mean(data))
    min_valley_depth <- mean(data)/10
    first_peak_found <- FALSE
    running_min <- data[left]
    candidate_right <- NULL
    candidate_score <- NULL

    right <- 2
    while (right < length(data)-1) {
        running_min <- min(running_min, data[right])
        is_local_peak <- (data[right - 1] < data[right]) && (data[right] >= data[right + 1])
        out_of_valley <- data[right] >= data[left]

        ## Out of valley check
        if(is_local_peak) {
            score <- min(data[left]-running_min, data[right]-running_min)
            message(paste("Considering valley", left, "-", right, "score:", score))
            if ((score > min_valley_depth)) {
                if (is.null(candidate_right)) {
                    message("Candidate valley!")
                    # valley is deep enough, so this is a valley
                    candidate_right <- right
                    candidate_score <- score
                } else if (data[candidate_right] < data[right]) {
                    message("Candidate valley!")
                    # valley is deep enough, so this is a valley
                    candidate_right <- right
                    candidate_score <- score
                } else {
                    message("passing on valley - candidate is better")
                    left <- right
                    running_min <- data[left]
                }
            } else if( (score <= min_valley_depth) && (data[left] < data[right]) ){
                message("Discarding valley")
                # valley is too shallow so we discard it
                left <- right
                running_min <- data[left]
                candidate_right <- NULL
                candidate_score <- NULL
            }

            if (data[left] < data[right]) {
                if (!is.null(candidate_right)) {
                    message("Accepting candidate ", left, "-", candidate_right)
                    result <- add_row(result, left = left, right = candidate_right, conv_depth = candidate_score)
                    left <- candidate_right
                    right <- candidate_right
                    running_min <- data[left]
                    candidate_right <- NULL
                    candidate_score <- NULL

                }
            }
        }

        if ((data[right-1] <= data[right]) && (data[right] > data[right+1])){
            # This discards the first 'valley' starting at position -2
            if (!first_peak_found) {
                message("Discarding prior to first peak")
                first_peak_found <- TRUE
                left <- right
                running_min <- data[left]
                candidate_right <- NULL
                candidate_score <- NULL
            }
        }
        right <- right + 1
    }

    if (!is.null(candidate_right)) {
        message("Accepting candidate ", left, "-", candidate_right)
        result <- add_row(result, left = left, right = candidate_right, conv_depth = score)
    }

    result <- result |>
        mutate(left = pmin(pmax(left - 2,1), len_gene), right = pmax(pmin(right - 2, len_gene),1))
    message("FINAL DF")
    message(result)

    return(result)
}

# Not used for valley score
#tp_peak_finder <- function(data) {
#  M <- max(data)
#  result <- tp_valley_finder(M - data)
#  return(result)
#}


#### Run smoothing convolution

cov_conv <- coverage |>
    group_by(gene_id) |>
    mutate(conv = get_conv_result(cov, conv_sd = CONV_SD))
cov_conv |> write_tsv(snakemake@output$smoothed_cov)

## Get all valleys (with Gauss smoothing)
valley_results <- list()
valley_results <- cov_conv |>
    group_by(across(all_of(group_vars))) |>
    group_modify(function(data, groups) { 
        message(paste("Looking at ", groups))
        alt_valley_finder(data$conv)
        #tp_valley_finder(data$conv)
        #tp_valley_finder_orig(data$conv)
    })

all_valley_results <- valley_results |>
  select(all_of(group_vars), left, right) |>
  mutate(optimum_type = "Valley")

#### Find local regions around Peaks and Valleys

# NOTE: this is not used in this script, which is just for computing the global valley score
#
### Helper function for finding intervals around valleys and peaks: from a vector of Boolean, extracts all the intervals of consecutive "TRUE"s
#get_interval_list <- function(mask){
#  interval_list <- list()
#  k <- 1
#  l <- 1
#  in_interval <- FALSE
#  if(mask[l] == TRUE){
#    in_interval <- TRUE
#  }else{
#    in_interval <- FALSE
#  }
#
#  mask <- c(mask, FALSE)
#
#  for(r in 1:length(mask)){
#    if(in_interval && (mask[r] == FALSE)){
#      interval_list[[k]] <- c(l,r-1)
#      k <- k+1
#      in_interval <- FALSE
#    }
#    else if(!in_interval && (mask[r] == TRUE)){
#      l <- r
#      in_interval <- TRUE
#    }
#  }
#  return(interval_list)
#}


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
local_valley_score |> write_tsv(snakemake@output$lvs)

transcript_valley_score <- local_valley_score |>
    left_join(coverage_stats, by=group_vars) |>
    group_by(across(all_of(group_vars))) |>
    summarize(
        transcript_valley_score = sum(lvs / mean_cov / transcript_length),
    )
transcript_valley_score |> write_tsv(snakemake@output$tvs)

if ("sample" %in% colnames(coverage)) {
    results <- transcript_valley_score |>
        group_by(sample) |>
        summarize(global_valley_score = sum(transcript_valley_score))
} else {
    results <- tibble(
        global_valley_score = sum(transcript_valley_score$transcript_valley_score),
    )
}
write_tsv(results, snakemake@output$gvs)
