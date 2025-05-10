
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
