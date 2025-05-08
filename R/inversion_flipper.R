
library(dplyr)
library(ggplot2)

read_anchorwave <- function(file_path) {
   #col_names <- c("rs", "re", "qs", "qe", "error", "qid", "rid", "strand")
   col_names <- c("rs", "re", "qs", "qe", "rid", "qid", "strand")
   data <- read.table(file_path, sep = "\t", header = TRUE, col.names = col_names)

   return(data)
}


find_matches <- function(data, min_size = 1000000) {
  # Calculate 'length' using 're' and 'rs'
  data <- as.data.frame(data)
  data$length <- data$re - data$rs

  # Create an empty list to store the results
  results <- list()

  # Loop through unique combinations of 'qid' and 'rid'
  unique_combinations <- unique(data[, c('qid', 'rid')])
  for (i in 1:nrow(unique_combinations)) {
    qid_i <- unique_combinations[i, 'qid']
    rid_i <- unique_combinations[i, 'rid']

    # Filter data for the current 'qid' and 'rid' combination
    filtered_data <- data[which(data$qid == qid_i & data$rid == rid_i), ]

    # Calculate the sum of 'length' for the current group
    total_length <- sum(filtered_data$length)

    # Number of concurrent alignments
    alignments <- nrow(filtered_data)

    # Create a result row
    result_row <- data.frame(qid = qid_i, rid = rid_i, total_length = total_length, alignments = alignments)

    # Add the result row to the results list
    results[[i]] <- result_row
  }

  # Combine the results into a data frame
  result_df <- do.call(rbind, results)

  # Order the result data frame by 'total_length' in descending order
  ordered_result <- result_df[order(result_df$total_length, decreasing = TRUE), ]
  ordered_result <- ordered_result[which(ordered_result$total_length >= min_size),]

  #Remove rows where there is only a single alignment for the comparison (therefore no inversions)
  ordered_result <- ordered_result[which(ordered_result$alignments > 1),]
  return(ordered_result)
}

rate_order_1 <- function(dataframe){
  order_1 <- dataframe[order(dataframe$rs),]
  problems <- 0
  for (i in 1:(nrow(order_1)-1)){
    problems <- problems+(order_1$qs[i+1] - order_1$qe[i])
  }
  return(problems)
}

rate_order_2 <- function(dataframe){
  order_1 <- dataframe[order(dataframe$rs),]
  problems <- 0
  gaps <- 0
  for (i in 1:(nrow(order_1)-1)){
    gaps <- gaps+abs((order_1$qs[i+1] - order_1$qe[i]))
    if (order_1$qs[i] > order_1$qs[i+1]){
      problems <- problems+1
    }
  }
  for (i in 1:(nrow(order_1))){
    if (order_1$qs[i] > order_1$qe[i]){
      problems <- problems+1
    }
  }
  return(c(problems, gaps))
}

check_overlap <- function(dataframe){
  order_1 <- dataframe[order(dataframe$qs),]
  problems <- 0
  for (i in 1:(nrow(order_1))){
    if (order_1$qs[i] > order_1$qe[i]){
      start <- order_1$qe[i]
      end <- order_1$qs[i]
      order_1$qe[i] <- end
      order_1$qs[i] <- start
    }
  }
  for (i in 1:(nrow(order_1)-1)){
    if (order_1$qe[i] > order_1$qs[i+1]){
      problems <- problems+1
    }
    if(order_1$qe[i] < 0){
      problems <- problems+1
    }
    if(order_1$qs[i] < 0){
      problems <- problems+1
    }
  }
  return(problems)
}

rate_direction <- function(dataframe){
  order_1 <- dataframe[order(dataframe$rs),]
  problems <- 0
  for (i in 1:(nrow(order_1)-1)){
    if (order_1$qs[i] > order_1$qe[i]){
      problems <- problems+1
    }
  }
  return(problems)
}


flip_orientation <- function(dataframe, boundary_start, boundary_end){
  if (nargs() <= 2){return()}
  all_edges <- c(dataframe$qs, dataframe$qe)
  if (boundary_start %in% all_edges & boundary_end %in% all_edges){

  }else{
    return()
  }
  mid <- (boundary_start + boundary_end) /2
  to_be_flipped <- dataframe[dataframe$qs >= boundary_start & dataframe$qe <= boundary_end ,]
  if (nrow(to_be_flipped) > 0){
    to_be_flipped$flipped <- T
  }

  not_flipped <- dataframe[!dataframe$rs %in% to_be_flipped$rs,]
  if (nrow(not_flipped) > 0){
    not_flipped$flipped <- F
  }
  flipped_rows <- to_be_flipped
  flipped_rows$qs <- mid - (to_be_flipped$qs -mid)
  flipped_rows$qe <- mid - (to_be_flipped$qe -mid)

  return_data <- rbind(not_flipped, flipped_rows)
  return(return_data)
}



full_flip_assess <- function(dataframe){
  original_score <- rate_order_2(dataframe)[1]

  all_breakpoints <- c(dataframe$qs, dataframe$qe)
  boundary_1 <- min(all_breakpoints)
  boundary_2 <- max(all_breakpoints)
  full_flip_score <- rate_order_2(flip_orientation(dataframe, boundary_1, boundary_2))[1]
  relative_score <- original_score - full_flip_score
  return(relative_score)
}

full_flip <- function(dataframe){

  all_breakpoints <- c(dataframe$qs, dataframe$qe)
  boundary_1 <- min(all_breakpoints)
  boundary_2 <- max(all_breakpoints)
  full_flipped <- flip_orientation(dataframe, boundary_1, boundary_2)
  return(full_flipped)
}


find_all_boundaries <- function(dataframe){
  tmp_dataframe <- dataframe
  boundaries_start <- NULL
  boundaries_end <- NULL
  for (i in 1:(nrow(tmp_dataframe))){
    if (tmp_dataframe$qs[i] < tmp_dataframe$qe[i]){
      boundaries_start <- c(boundaries_start, tmp_dataframe$qs[i])
      boundaries_end <- c(boundaries_end, tmp_dataframe$qe[i])
    }else{
      boundaries_start <- c(boundaries_start, tmp_dataframe$qe[i])
      boundaries_end <- c(boundaries_end, tmp_dataframe$qs[i])
    }
  }
  good_flips <- data.frame()
  n<-1
  for (i in 1:length(boundaries_start)){
    for (j in 1:length(boundaries_end)){
      boundary_1 <- boundaries_start[i]
      boundary_2 <- boundaries_end[j]
      if(boundary_1 >= boundary_2){next}
      overlap_score <- check_overlap(flip_orientation(dataframe, boundary_1, boundary_2))
      if (overlap_score > 0){
        next
      }
      tmp <- data.frame(boundary_1 = boundary_1, boundary_2 = boundary_2, n=n)
      n<- n+1
      good_flips <- rbind(good_flips, tmp)
    }
  }
  return(good_flips)
}



find_all_boundaries_deep<-function(dataframe){
  boundaries <- data.frame()
  level_1 <- find_all_boundaries(dataframe)
  for (i in 1:nrow(level_1)){
    print(paste0("Level 1 ",i))
    boundary_1.1 <- level_1$boundary_1[i]
    boundary_2.1 <- level_1$boundary_2[i]
    level_1_data <- flip_orientation(dataframe, boundary_1.1, boundary_2.1)
    save_results <- data.frame(boundary_1 = c(boundary_1.1),
                               boundary_2 = c(boundary_2.1),
                               level = 1, id= paste0(i))
    boundaries <- rbind(boundaries,save_results)
    level_2 <- find_all_boundaries(level_1_data)
    for (j in 1:nrow(level_2)){

      boundary_1.2 <- level_2$boundary_1[j]
      boundary_2.2 <- level_2$boundary_2[j]
      level_2_data <- flip_orientation(level_1_data, boundary_1.2, boundary_2.2)
      save_results <- data.frame(boundary_1 = c(boundary_1.1, boundary_1.2),
                                 boundary_2 = c(boundary_2.1, boundary_2.2),
                                 level = 2, id= paste0(i,"-",j))
      boundaries <- rbind(boundaries,save_results)
      level_3 <- find_all_boundaries(level_2_data)
      for (k in 1:nrow(level_3)){

        #Save all flips
        boundary_1.3 <- level_3$boundary_1[k]
        boundary_2.3 <- level_3$boundary_2[k]
        save_results <- data.frame(boundary_1 = c(boundary_1.1, boundary_1.2, boundary_1.3),
                                   boundary_2 = c(boundary_2.1, boundary_2.2, boundary_2.3),
                                   level = 3, id= paste0(i,"-",j,"-",k))
        boundaries <- rbind(boundaries,save_results)

      }
    }
  }
  return(boundaries)
}

level_1_test <- function(dataframe){
  level_1 <- find_all_boundaries(dataframe)
  boundaries <- data.frame()
  for (i in 1:nrow(level_1)){
    #print(paste0("Level 1 ",i))
    boundary_1.1 <- level_1$boundary_1[i]
    boundary_2.1 <- level_1$boundary_2[i]
    level_1_flip <- flip_orientation(dataframe, boundary_1.1, boundary_2.1)

    overlap_score <- check_overlap(level_1_flip)
    new_stretches <- test_for_stretches(level_1_flip)
    if (overlap_score > 0){
      next
    }

    save_results <- data.frame(boundary_1 = c(boundary_1.1),
                               boundary_2 = c(boundary_2.1),
                               level = 1, id= paste0(i),
                               stretches=new_stretches)
    boundaries <- rbind(boundaries,save_results)

  }
  return(boundaries)
}

level_2_test <- function(dataframe){
  level_1 <- find_all_boundaries(dataframe)
  boundaries <- data.frame()
  for (i in 1:nrow(level_1)){
    #print(paste0("Level 1 ",i))
    boundary_1.1 <- level_1$boundary_1[i]
    boundary_2.1 <- level_1$boundary_2[i]
    level_1_flip <- flip_orientation(dataframe, boundary_1.1, boundary_2.1)

    overlap_score <- check_overlap(level_1_flip)
    new_stretches <- test_for_stretches(level_1_flip)
    if (overlap_score > 0){
      next
    }

    save_results <- data.frame(boundary_1 = c(boundary_1.1),
                               boundary_2 = c(boundary_2.1),
                               level = 1, id= paste0(i),
                               stretches=new_stretches)
    boundaries <- rbind(boundaries,save_results)
    level_2 <- level_1[-i,]
    if (nrow(level_2) == 0){next;}
    for (j in 1:nrow(level_2)){

      boundary_1.2 <- level_2$boundary_1[j]
      boundary_2.2 <- level_2$boundary_2[j]
      level_2_flip <- flip_orientation(level_1_flip, boundary_1.2, boundary_2.2)
      if (is.null(level_2_flip)){next}
      overlap_score <- check_overlap(level_2_flip)
      new_stretches <- test_for_stretches(level_2_flip)
      if (overlap_score > 0){
        next
      }


      save_results <- data.frame(boundary_1 = c(boundary_1.1, boundary_1.2),
                                 boundary_2 = c(boundary_2.1, boundary_2.2),
                                 level = 2, id= paste0(i,"-",j),
                                 stretches=new_stretches)
      boundaries <- rbind(boundaries,save_results)

    }
  }
  return(boundaries)
}

level_3_test <- function(dataframe){
  level_1 <- find_all_boundaries(dataframe)
  boundaries <- data.frame()
  for (i in 1:nrow(level_1)){
    #print(paste0("Level 1 ",i))
    boundary_1.1 <- level_1$boundary_1[i]
    boundary_2.1 <- level_1$boundary_2[i]
    level_1_flip <- flip_orientation(dataframe, boundary_1.1, boundary_2.1)

    overlap_score <- check_overlap(level_1_flip)
    new_stretches <- test_for_stretches(level_1_flip)
    if (overlap_score > 0){
      next
    }

    save_results <- data.frame(boundary_1 = c(boundary_1.1),
                               boundary_2 = c(boundary_2.1),
                               level = 1, id= paste0(i),
                               stretches=new_stretches)
    boundaries <- rbind(boundaries,save_results)
    level_2 <- level_1[-i,]
    if (nrow(level_2) == 0){next;}
    for (j in 1:nrow(level_2)){

      boundary_1.2 <- level_2$boundary_1[j]
      boundary_2.2 <- level_2$boundary_2[j]
      level_2_flip <- flip_orientation(level_1_flip, boundary_1.2, boundary_2.2)
      if (is.null(level_2_flip)){next}
      overlap_score <- check_overlap(level_2_flip)
      new_stretches <- test_for_stretches(level_2_flip)
      if (overlap_score > 0){
        next
      }


      save_results <- data.frame(boundary_1 = c(boundary_1.1, boundary_1.2),
                                 boundary_2 = c(boundary_2.1, boundary_2.2),
                                 level = 2, id= paste0(i,"-",j),
                                 stretches=new_stretches)
      boundaries <- rbind(boundaries,save_results)
      level_3 <- level_2[-j,]
      if (nrow(level_3) == 0){next;}
      for (k in 1:nrow(level_3)){

        #Save all flips
        boundary_1.3 <- level_3$boundary_1[k]
        boundary_2.3 <- level_3$boundary_2[k]
        level_3_flip <- flip_orientation(level_2_flip, boundary_1.3, boundary_2.3)
        if (is.null(level_3_flip)){next}

        overlap_score <- check_overlap(level_3_flip)
        new_stretches <- test_for_stretches(level_3_flip)
        if (overlap_score > 0){
          next
        }
        save_results <- data.frame(boundary_1 = c(boundary_1.1, boundary_1.2, boundary_1.3),
                                   boundary_2 = c(boundary_2.1, boundary_2.2, boundary_2.3),
                                   level = 3, id= paste0(i,"-",j,"-",k),
                                   stretches=new_stretches)
        boundaries <- rbind(boundaries,save_results)

      }
    }
  }
  return(boundaries)
}




test_all_flips <- function(dataframe){
  original_dataframe <- dataframe



  #Remove initial overlaps
  tmp_dataframe <- remove_overlap(dataframe)
  if (!is.null(tmp_dataframe)){
    dataframe <- tmp_dataframe
  }
  #Do initial merging
  merging_ongoing <- T
  while (merging_ongoing){
    if(nrow(dataframe) == 1){break}
    tmp_dataframe <- find_stretches(dataframe)
    if (!is.null(tmp_dataframe)){
      dataframe <- tmp_dataframe
      #print("merged stretch")
    }else{
      merging_ongoing <- F
    }
  }
  if(nrow(dataframe) == 1){return()}

  round<-1
  inversion_list <- data.frame()
  while(nrow(dataframe) > 1){
    round <- round+1
    boundaries <- level_1_test(dataframe)
    #print(paste0("Round ", round, ": Trying level 1"))
    if (max(boundaries$stretches) == 0){
      #print(paste0("Round ", round, ": Trying level 2"))
      boundaries <- level_2_test(dataframe)
      if (max(boundaries$stretches) == 0){
        #print(paste0("Round ", round, ": Trying level 3"))
        boundaries <- level_3_test(dataframe)
      }
    }

    #pick set with lowest number of problems and fewest flips, then randomly select
    top_choices <- boundaries[which(boundaries$stretches == max(boundaries$stretches)),]
    top_choices <- top_choices[which(top_choices$level == min(top_choices$level)),]
    chosen_flip <- top_choices[sample(1:nrow(top_choices)), ][1,]

    n_flips <- boundaries[which(boundaries$id == chosen_flip$id[1]),]$level[1]
    flips_testing <- boundaries[which(boundaries$id == chosen_flip$id[1]),]
    for (i in 1:n_flips){
      track_inversion <- data.frame(boundary_1= flips_testing[i,1], boundary_2=flips_testing[i,2])
      inversion_list <- rbind(inversion_list, track_inversion)
      dataframe <- flip_orientation(dataframe, flips_testing[i,1], flips_testing[i,2])
    }
    #Do a round of merging
    merging_ongoing <- T
    while (merging_ongoing){
      if(nrow(dataframe) == 1){break}
      tmp_dataframe <- find_stretches(dataframe)
      if (!is.null(tmp_dataframe)){
        dataframe <- tmp_dataframe
        #print("merged stretch")
      }else{
        merging_ongoing <- F
      }
    }
  }
  return(inversion_list)
}




find_stretches <- function(dataframe){
  order_1 <- dataframe[order(dataframe$rs),]
  for (i in 1:(nrow(order_1)-1)){
    #Check if they're both in the right orientation
    if (order_1$qs[i] > order_1$qe[i]){
      next
    }
    if (order_1$qs[i+1] > order_1$qe[i+1]){
      next
    }
    #Check if they're in the right order.
    if (order_1$qe[i] > order_1$qe[i+1]){
      next
    }

    #Check if they are consecutive
    gap_start <- order_1$qe[i]
    gap_end <- order_1$qs[i+1]
    if (sum(order_1$qs > gap_start & order_1$qs < gap_end) > 0){
      next
    }
    new_row <- data.frame()
    #Merge the consecutive windows.
    if(ncol(dataframe) == 5){
      new_row <- data.frame(rs = order_1$rs[i],
                            re = order_1$re[i+1],
                            qs = order_1$qs[i],
                            qe = order_1$qe[i+1],
                            flipped=NA)
    }else{
      new_row <- data.frame(rs = order_1$rs[i],
                            re = order_1$re[i+1],
                            qs = order_1$qs[i],
                            qe = order_1$qe[i+1])
    }

    dataframe <- order_1[-c(i,i+1),]
    dataframe <- rbind(new_row, dataframe)
    return(dataframe)
  }
  #Looking for reverse stretches and merge those
  for (i in 1:(nrow(order_1)-1)){
    #Check if they're both in the reverse orientation
    if (order_1$qs[i] < order_1$qe[i]){
      next
    }
    if (order_1$qs[i+1] < order_1$qe[i+1]){
      next
    }
    #Check if they're in the right order.
    if (order_1$qe[i] < order_1$qe[i+1]){
      next
    }

    #Check if they are consecutive
    gap_start <- order_1$qe[i]
    gap_end <- order_1$qs[i+1]
    if (sum(order_1$qs < gap_start & order_1$qs > gap_end) > 0){
      next
    }
    new_row <- data.frame()
    #Merge the consecutive windows.
    if(ncol(dataframe) == 5){
      new_row <- data.frame(rs = order_1$rs[i],
                            re = order_1$re[i+1],
                            qs = order_1$qs[i],
                            qe = order_1$qe[i+1],
                            flipped=NA)
    }else{
      new_row <- data.frame(rs = order_1$rs[i],
                            re = order_1$re[i+1],
                            qs = order_1$qs[i],
                            qe = order_1$qe[i+1])
    }

    dataframe <- order_1[-c(i,i+1),]
    dataframe <- rbind(new_row, dataframe)
    return(dataframe)
  }
}

test_for_stretches <- function(dataframe){
  order_1 <- dataframe[order(dataframe$rs),]
  stretches <- 0
  for (i in 1:(nrow(order_1)-1)){
    #Check if they're both in the right orientaion
    if (order_1$qs[i] > order_1$qe[i]){
      next
    }
    if (order_1$qs[i+1] > order_1$qe[i+1]){
      next
    }
    #Check if they're in the right order.
    if (order_1$qe[i] > order_1$qe[i+1]){
      next
    }
    #Check if they are consecutive
    gap_start <- order_1$qe[i]
    gap_end <- order_1$qs[i+1]
    if (sum(order_1$qs > gap_start & order_1$qs < gap_end) > 0){
      next
    }
    #Merge the consecutive windows.
    stretches <- stretches + 1
  }
  return(stretches)
}

remove_overlap <- function(dataframe) {
  order_1 <- dataframe[order(dataframe$rs),]
  rows_to_remove <- c()

  repeat {
    overlaps_found <- FALSE

    for (i in 1:(nrow(order_1)-1)) {
      if (order_1$re[i] > order_1$rs[i+1]) {
        size_1 <- order_1$re[i] - order_1$rs[i]
        size_2 <- order_1$re[i+1] - order_1$rs[i+1]

        if (size_1 > size_2) {
          rows_to_remove <- c(rows_to_remove, i + 1)
        } else {
          rows_to_remove <- c(rows_to_remove, i)
        }

        overlaps_found <- TRUE
      }
    }

    if (overlaps_found) {
      order_1 <- order_1[-rows_to_remove, ]
      rows_to_remove <- c()
    } else {
      break
    }
  }

  cleaned_dataframe <- order_1
  return(cleaned_dataframe)
}



create_flipped_alignment <- function(alignment, flips){
  tmp_dataframe <- alignment
  for (i in 1:nrow(flips)){
    tmp_dataframe <- flip_orientation(tmp_dataframe, flips[i,3], flips[i,4])
  }
  return(tmp_dataframe)
}




plot_flip_schedule <- function(alignment,flips){
  all_alignments <- data.frame()
  original <- alignment
  original$flips <- 0
  original$flipped <- F
  all_alignments <- rbind(all_alignments, original)
  for (i in 1:nrow(flips)){
    #print(i)
    alignment <- flip_orientation(alignment, flips[i,1],flips[i,2])
    flipped <- alignment
    flipped$flips <- i
    all_alignments <- rbind(all_alignments, flipped)


  }
  return(all_alignments)
}

affected_windows <- function(dataframe, boundary_start, boundary_end){

  to_be_flipped <- dataframe[dataframe$qs >= boundary_start & dataframe$qe <= boundary_end ,]
  to_be_flipped <- to_be_flipped[,c(1,2)]
  return(to_be_flipped)
}

list_flipped_windows <-  function(alignment, flips){
  tmp_dataframe <- alignment
  flipped_windows <- data.frame()
  for (i in 1:nrow(flips)){
    tmp <- affected_windows(tmp_dataframe, flips[i,1], flips[i,2])
    tmp$inversion_n <- i
    flipped_windows <- rbind(flipped_windows, tmp)
    tmp_dataframe <- flip_orientation(tmp_dataframe, flips[i,1], flips[i,2])
  }
  return(flipped_windows)
}


basic_alignment_plot <- function(data, x_col1 = "rs", x_col2 = "re",
                                 y_col1 = "qs", y_col2 = "qe", xlab = "Reference (Mbp)", ylab = "Query (Mbp)", main = NULL) {
  # Create a new plot
  plot(NA, NA, xlim = range(data[[x_col1]], data[[x_col2]]) / 1000000,
       ylim = range(data[[y_col1]], data[[y_col2]]) / 1000000,
       xlab = xlab, ylab = ylab, main = main)

  # Add segments to the plot
  segments(x0 = data[[x_col1]] / 1000000, x1 = data[[x_col2]] / 1000000,
           y0 = data[[y_col1]] / 1000000, y1 = data[[y_col2]] / 1000000)
}

multi_iterate <- function(input_data, iterations=100, return="best_one",verbose=F){
  all_attempts <- data.frame()

  for (i in 1:iterations){
    if(verbose){print(paste("Trying iteration",i))}
    result <- test_all_flips(input_data)
    if (is.null(result)){next}
    result$replicate <- i
    result$n_flips <- nrow(result)
    all_attempts <- rbind(all_attempts, result)
  }
  if(nrow(all_attempts) == 0){return()}
  top_attempts <- all_attempts[which(all_attempts$n_flips == min(all_attempts$n_flips)),]
  top_attempt <- top_attempts[which(top_attempts$replicate == min(top_attempts$replicate)),]
  if(return == "all"){
    return(all_attempts)
  }else if(return == "best_all"){
    return(top_attempts)
  }else if (return == "best_one"){
    return(top_attempt)
  }else{
    warning("Return needs to be 'all', 'best_all', or 'best_one'")
  }
}


plot_final_series <- function(data, flips, ...) {
  flips_plotting <- flips
  flips_plotting$flips <- 1:nrow(flips_plotting)
  min_x = min(data$rs)

  # Create a list of additional arguments to pass to ggplot
  ggplot_args <- list(...)

  # Start building the ggplot
  p <- ggplot(data, ) +
    geom_segment(ggplot2::aes(x = rs/1000000, xend = re/1000000, y = qs/1000000, yend = qe/1000000, color = flipped)) +
    xlab("Reference (Mbp)") + ylab("Query (Mbp)") +
    scale_color_manual(values = c("dark blue", "red"), name = "Flipped",
                       labels=c("No","Yes")) +
    geom_segment(data = flips_plotting, aes(y = boundary_1/1000000, yend = boundary_2/1000000, x = min_x/1000000, xend = min_x/1000000)) +
    facet_wrap(~flips)

  # Add any additional arguments to ggplot
  if (!missing(...)) {
    p <- p + ...
  }
  return(p)
}


find_all_inversions <- function(file_name, prefix="deverter.out", min_size=10000000,iterations=100){
  print(paste("Reading", file_name))
  data <- read_anchorwave(file_name)
  print(paste("Finding all chromosome-chromosome matches with a minimum size of",min_size))
  matches <- find_matches(data, min_size)
  print(paste("Found",nrow(matches), "chromosome matches"))
  inversion_regions = data.frame()
  inversion_boundaries = data.frame()
  pdf(paste0(prefix,".pdf"))
  for (i in 1:nrow(matches)){
    rmatch <- matches[i,2]
    qmatch <- matches[i,1]
    print(paste0("Match #", i," Comparing ",rmatch, " to ",qmatch))
    test <- data[which(data$qid == qmatch & data$rid == rmatch),c(1:4)]
    result <- multi_iterate(test, iterations=iterations)
    if(is.null(result)){next}
    flipped_windows <- list_flipped_windows(test, result)
    flipped_windows$qid <- qmatch
    flipped_windows$rid <- rmatch
    flipped_windows$inversion_id <- paste(flipped_windows$rid, flipped_windows$qid,"Inversion", flipped_windows$inversion_n, sep = ".")
    inversion_regions <- rbind(inversion_regions, flipped_windows)
    plotting_result <- plot_flip_schedule(test, result)

    saved_result <- result[,3:4]
    saved_result$qid <- qmatch
    saved_result$rid <- rmatch
    inversion_boundaries <- rbind(inversion_boundaries, saved_result)
    plot(
      plot_final_series(plotting_result, result,ggtitle(paste0(rmatch, "-", qmatch)))
    )
  }
  dev.off()
  write.table(inversion_regions, paste0(prefix,".txt"),quote=F, row.names = F)
}


simulate_inversions <- function(size=10000,inversions=1){

  original <- data.frame()
  #Variables just controlling the little windows used to find inversions.
  window_size <- 20
  window_length <- window_size/2
  for (i in seq(window_size, size-window_size, window_size)){
    row <- data.frame(rs=i,re=i+window_length,qs=i,qe=i+window_length)
    original <- rbind(original, row)
  }
  flipped <- original
  breakpoints <- data.frame()
  for (i in 1:inversions ){
    b1 <- 2
    b2 <- 1
    while(b1 > b2){
      b1 <- sample(original$rs, 1)
      b2 <- sample(original$re, 1)
    }
    tmp <- data.frame(start=b1,end=b2)
    breakpoints <- rbind(breakpoints, tmp)
  }

  #Genome is broken up enough, now to rotate it
  for (i in 1:inversions){
    flipped <- flip_orientation(flipped, breakpoints$start[i], breakpoints$end[i])
  }
  output <- list(breakpoints,flipped)
  return(output)
}

benchmark_sim <- function(inversions=1, iterations=10,replicates=1){
  benchmark_result <- data.frame()
  for (x in 1:replicates){
    test_sim <- simulate_inversions(inversions=inversions)
    estimate <- multi_iterate(test_sim[[2]],iterations=iterations)
    smushed_real <- paste(test_sim[[1]]$start, test_sim[[1]]$end)
    smushed_estimate <-  paste(estimate$boundary_1, estimate$boundary_2)
    #Check if they're right.
    #Check if real inversions were found.
    inversions_found <- 0
    for (i in 1:length(smushed_real)){
      if (smushed_real[i] %in% smushed_estimate ){
        inversions_found <- inversions_found + 1
      }
    }
    missed_inversions <- inversions - inversions_found
    #Check if extra inversions were found
    estimates_found <- 0
    for (i in 1:length(smushed_estimate)){
      if (smushed_estimate[i] %in% smushed_real ){
        estimates_found <- estimates_found + 1
      }
    }
    extra_inversions <- length(smushed_estimate) - estimates_found
    result <- data.frame(rep=x,
                         missed_inversions=missed_inversions,
                         extra_inversions=extra_inversions,
                         total_inversions=inversions,
                         iterations=iterations)
    benchmark_result <- rbind(benchmark_result, result)
  }
  return(benchmark_result)

}


