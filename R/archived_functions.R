test_all_flips <- function(dataframe){
  base_score <- rate_order(dataframe)
  tmp_dataframe <- dataframe
  for (i in 1:(nrow(tmp_dataframe))){
    if (tmp_dataframe$qs[i] > tmp_dataframe$qe[i]){
      start <- tmp_dataframe$qe[i]
      end <- tmp_dataframe$qs[i]
      tmp_dataframe$qe[i] <- end
      tmp_dataframe$qs[i] <- start
    }
  }
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
  for (i in 1:length(boundaries_start)){
    for (j in 1:length(boundaries_end)){
      boundary_1 <- boundaries_start[i]
      boundary_2 <- boundaries_end[j]
      if(boundary_1 >= boundary_2){next}
      tmp_score <- rate_order(flip_orientation(dataframe, boundary_1, boundary_2))  - base_score
      overlap_score <- check_overlap(flip_orientation(dataframe, boundary_1, boundary_2))
      new_stretches <- test_for_stretches(flip_orientation(dataframe, boundary_1, boundary_2))
      if (overlap_score > 0){
        next
      }
      tmp <- data.frame(new_stretches = new_stretches, score=tmp_score, boundary_1 = boundary_1, boundary_2 = boundary_2)
      good_flips <- rbind(good_flips, tmp)
    }
  }
  return(good_flips)
}

test_all_double_flips <- function(dataframe){
  base_score <- rate_order(dataframe)
  tmp_dataframe <- dataframe
  for (i in 1:(nrow(tmp_dataframe))){
    if (tmp_dataframe$qs[i] > tmp_dataframe$qe[i]){
      start <- tmp_dataframe$qe[i]
      end <- tmp_dataframe$qs[i]
      tmp_dataframe$qe[i] <- end
      tmp_dataframe$qs[i] <- start
    }
  }

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
  for (i in 1:length(boundaries_start)){
    for (j in 1:length(boundaries_end)){
      for (i2 in 1:length(boundaries_start)){
        for (j2 in 1:length(boundaries_end)){
          boundary_1 <- boundaries_start[i]
          boundary_2 <- boundaries_end[j]
          boundary_3 <- boundaries_start[i2]
          boundary_4 <- boundaries_end[j2]
          if(boundary_1 == boundary_3 & boundary_2 == boundary_4){next}
          if(boundary_1 >= boundary_2){next}
          if(boundary_3 >= boundary_4){next}
          #print(paste("Trying", boundary_1, boundary_2, boundary_3, boundary_4))
          flipped_data_frame <- flip_orientation(
            flip_orientation(
              dataframe, boundary_1, boundary_2
            ), boundary_3, boundary_4)
          tmp_score <- rate_order(flipped_data_frame)  - base_score
          overlap_score <- check_overlap(flipped_data_frame)
          new_stretches <- test_for_stretches(flipped_data_frame)
          if (overlap_score > 0){
            next
          }
          if (new_stretches == 0){
            next
          }
          tmp <- data.frame(new_stretches = new_stretches, score=tmp_score, boundary_1 = boundary_1, boundary_2 = boundary_2,
                            boundary_3 = boundary_3, boundary_4 = boundary_4)
          good_flips <- rbind(good_flips, tmp)
        }
      }
    }
  }
  return(good_flips)
}

test_all_triple_flips <- function(dataframe){
  base_score <- rate_order(dataframe)
  tmp_dataframe <- dataframe
  for (i in 1:(nrow(tmp_dataframe))){
    if (tmp_dataframe$qs[i] > tmp_dataframe$qe[i]){
      start <- tmp_dataframe$qe[i]
      end <- tmp_dataframe$qs[i]
      tmp_dataframe$qe[i] <- end
      tmp_dataframe$qs[i] <- start
    }
  }

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
  for (i1 in 1:length(boundaries_start)){
    for (j1 in 1:length(boundaries_end)){
      for (i2 in 1:length(boundaries_start)){
        for (j2 in 1:length(boundaries_end)){
          for (i3 in 1:length(boundaries_start)){
            for (j3 in 1:length(boundaries_end)){
              boundary_1 <- boundaries_start[i1]
              boundary_2 <- boundaries_end[j1]
              boundary_3 <- boundaries_start[i2]
              boundary_4 <- boundaries_end[j2]
              boundary_5 <- boundaries_start[i3]
              boundary_6 <- boundaries_end[j3]
              if(boundary_1 >= boundary_2){next}
              if(boundary_3 >= boundary_4){next}
              if(boundary_5 >= boundary_6){next}
              flipped_data_frame <- flip_orientation(
                flip_orientation(
                  flip_orientation(dataframe,
                                   boundary_1, boundary_2),
                  boundary_3, boundary_4),
                boundary_5, boundary_6)
              tmp_score <- rate_order(flipped_data_frame)  - base_score
              overlap_score <- check_overlap(flipped_data_frame)
              new_stretches <- test_for_stretches(flipped_data_frame)
              if (overlap_score > 0){
                next
              }
              if (new_stretches == 0){
                next
              }
              tmp <- data.frame(new_stretches = new_stretches, score=tmp_score, boundary_1 = boundary_1, boundary_2 = boundary_2,
                                boundary_3 = boundary_3, boundary_4 = boundary_4,boundary_5,boundary_6)
              good_flips <- rbind(good_flips, tmp)
            }
          }
        }
      }
    }
  }
  return(good_flips)
}


iterate_structural_flips <- function(dataframe, verbose=F){
  #Remove initial overlaps
  tmp_dataframe <- remove_overlap(dataframe)
  if (!is.null(tmp_dataframe)){
    dataframe <- tmp_dataframe
  }
  #Do initial merging
  merging_ongoing <- T
  while (merging_ongoing){
    if(nrow(dataframe == 1)){break}
    tmp_dataframe <- find_stretches(dataframe)
    if (!is.null(tmp_dataframe)){
      dataframe <- tmp_dataframe
      print("merged stretch")
    }else{
      merging_ongoing <- F
    }
  }
  if(nrow(dataframe) == 1){return()}
  improved_possible <- 1
  loops <- 0
  flips_done <- data.frame()
  #Try single flips
  while (improved_possible >= 1){
    if (nrow(dataframe) == 1){break}
    loops <- loops +1
    if(verbose == T){print(paste("Trying single flips", loops))}
    better_scores <- test_all_flips(dataframe)
    try_double <- F
    try_triple <- F
    if (nrow(better_scores) == 0 ){

    }else{
      top_scores <- max(better_scores$new_stretches)
      if (top_scores == 0){
        try_double <- T
      }
    }


    if (try_double){
      if(verbose == T){print(paste("Trying double flips", loops))}
      better_scores <- test_all_double_flips(dataframe)
      try_triple <- F
      if (nrow(better_scores) == 0 ){
        try_triple <- T
      }else{
        top_scores <- max(better_scores$new_stretches)
        if (top_scores == 0){
          try_triple <- T
        }
      }
      if (try_triple){
        if(verbose == T){print(paste("Trying triple flips", loops))}
        #Try three flips
        better_scores <- test_all_triple_flips(dataframe)
        if (nrow(better_scores) == 0){break}
        top_scores <- max(better_scores$new_stretches)
        if (top_scores == 0){break}
        better_scores <- better_scores[which(better_scores$new_stretches == top_scores),]
        rows <- sample(nrow(better_scores))
        flip_chosen <- head(better_scores[rows, ],1)
        dataframe <- flip_orientation(
          flip_orientation(
            flip_orientation(
              dataframe, flip_chosen[,3],flip_chosen[,4]),
            flip_chosen[,5],flip_chosen[,6]),
          flip_chosen[,7],flip_chosen[,8])
        flip_chosen <- data.frame(new_stretches = c(head(better_scores[rows, ],1)[,1],
                                                    head(better_scores[rows, ],1)[,1],
                                                    head(better_scores[rows, ],1)[,1]),
                                  score =  c(head(better_scores[rows, ],1)[,2],
                                             head(better_scores[rows, ],1)[,2],
                                             head(better_scores[rows, ],1)[,2]),
                                  boundary_1 = c(head(better_scores[rows, ],1)[,3],
                                                 head(better_scores[rows, ],1)[,5],
                                                 head(better_scores[rows, ],1)[,7]),
                                  boundary_2 = c(head(better_scores[rows, ],1)[,4],
                                                 head(better_scores[rows, ],1)[,6],
                                                 head(better_scores[rows, ],1)[,8]))
        flips_done <- rbind(flips_done, flip_chosen)
        if (nrow(dataframe) == 1){break}
        merging_ongoing <- T
        while (merging_ongoing){
          if (nrow(dataframe) == 1){break}
          tmp_dataframe <- find_stretches(dataframe)
          if (!is.null(tmp_dataframe)){
            dataframe <- tmp_dataframe
            if(verbose == T){print("merged stretch")}
          }else{
            merging_ongoing <- F
          }
        }
      }else{
        top_scores <- max(better_scores$new_stretches)
        if (top_scores == 0){break}
        better_scores <- better_scores[which(better_scores$new_stretches == top_scores),]
        rows <- sample(nrow(better_scores))
        flip_chosen <- head(better_scores[rows, ],1)
        dataframe <- flip_orientation(flip_orientation(dataframe, flip_chosen[,3],flip_chosen[,4]), flip_chosen[,5],flip_chosen[,6])
        flip_chosen <- data.frame(new_stretches = c(head(better_scores[rows, ],1)[,1], head(better_scores[rows, ],1)[,1]),
                                  score =  c(head(better_scores[rows, ],1)[,2], head(better_scores[rows, ],1)[,2]),
                                  boundary_1 = c(head(better_scores[rows, ],1)[,3], head(better_scores[rows, ],1)[,5]),
                                  boundary_2 = c(head(better_scores[rows, ],1)[,4], head(better_scores[rows, ],1)[,6]))
        flips_done <- rbind(flips_done, flip_chosen)
        if (nrow(dataframe) == 1){break}
        merging_ongoing <- T
        while (merging_ongoing){
          if (nrow(dataframe) == 1){break}
          tmp_dataframe <- find_stretches(dataframe)
          if (!is.null(tmp_dataframe)){
            dataframe <- tmp_dataframe
            if(verbose == T){print("merged stretch")}
          }else{
            merging_ongoing <- F
          }
        }
      }
    }else{

      better_scores <- better_scores[which(better_scores$new_stretches == top_scores),]
      rows <- sample(nrow(better_scores))
      flip_chosen <- head(better_scores[rows, ],1)
      flips_done <- rbind(flips_done, flip_chosen)
      dataframe <- flip_orientation(dataframe, flip_chosen[,3],flip_chosen[,4])
      if (nrow(dataframe) == 1){break}
      merging_ongoing <- T
      while (merging_ongoing){
        if (nrow(dataframe) == 1){break}
        tmp_dataframe <- find_stretches(dataframe)
        if (!is.null(tmp_dataframe)){
          dataframe <- tmp_dataframe
          if(verbose == T){print("merged stretch")}
        }else{
          merging_ongoing <- F
        }
      }
    }


  }
  reported_dataframe <- dataframe
  reported_flips <- flips_done
  return(reported_flips)
}


test_all_singletons <- function(dataframe){
  base_score <- rate_direction(dataframe)
  tmp_dataframe <- dataframe
  for (i in 1:(nrow(tmp_dataframe))){
    if (tmp_dataframe$qs[i] > tmp_dataframe$qe[i]){
      start <- tmp_dataframe$qe[i]
      end <- tmp_dataframe$qs[i]
      tmp_dataframe$qe[i] <- end
      tmp_dataframe$qs[i] <- start
    }
  }
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
  for (i in 1:length(boundaries_start)){
    boundary_1 <- boundaries_start[i]
    boundary_2 <- boundaries_end[i]
    tmp_score <- rate_direction(flip_orientation(dataframe, boundary_1, boundary_2))  - base_score
    overlap_score <- check_overlap(flip_orientation(dataframe, boundary_1, boundary_2))
    if (overlap_score > 0){
      next
    }
    if (tmp_score < 0){
      tmp <- data.frame(score=tmp_score, boundary_1 = boundary_1, boundary_2 = boundary_2)
      good_flips <- rbind(good_flips, tmp)
    }
  }
  return(good_flips)
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







test_all_flips <- function(dataframe){
  original_dataframe <- dataframe

  full_flip_check <- full_flip_assess(dataframe)

  inversion_list <- data.frame()
  if (full_flip_check > 0){
    all_breakpoints <- c(dataframe$qs, dataframe$qe)
    boundary_1 <- min(all_breakpoints)
    boundary_2 <- max(all_breakpoints)
    full_flipped <- flip_orientation(dataframe, boundary_1, boundary_2)
    dataframe <- full_flipped
    track_inversion <- data.frame(boundary_1=boundary_1, boundary_2=boundary_2)
    inversion_list <- rbind(inversion_list, track_inversion)
  }


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


  round<-1
  while(nrow(dataframe) > 1){
    print(paste0("Round ",round))
    round <- round+1
    starting_score <- rate_order_2(dataframe)[1]
    possible_flips <- find_all_boundaries_deep(dataframe)
    possible_ids <- unique(possible_flips$id)
    tried_flips <- data.frame()
    for (x in possible_ids){
      n_flips <- possible_flips[which(possible_flips$id == x),]$level[1]
      flips_testing <- possible_flips[which(possible_flips$id == x),]
      tmp_dataframe <- dataframe
      for (i in 1:n_flips){
        tmp_dataframe <- flip_orientation(tmp_dataframe, flips_testing[i,1], flips_testing[i,2])
      }
      new_scores <- rate_order_2(tmp_dataframe)
      score_improvement <- starting_score - new_scores[1]
      tmp_score <- data.frame(id=x,score=score_improvement,flips_done=n_flips,gaps=new_scores[2])
      tried_flips <- rbind(tried_flips, tmp_score)
    }
    #pick set with lowest number of problems and fewest flips, then randomly select
    top_choices <- tried_flips[which(tried_flips$score == max(tried_flips$score)),]
    top_choices <- top_choices[which(top_choices$flips_done == min(top_choices$flips_done)),]
    top_choices <- top_choices[which(top_choices$gaps == min(top_choices$gaps)),]
    chosen_flip <- sample(top_choices)[1,]

    n_flips <- possible_flips[which(possible_flips$id == chosen_flip$id[1]),]$level[1]
    flips_testing <- possible_flips[which(possible_flips$id == chosen_flip$id[1]),]
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


