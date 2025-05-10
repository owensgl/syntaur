#Benchmarking this program
library(foreach)
library(doParallel)

# Set the number of cores to use (adjust as needed)
num_cores <- 6
cl <- makeCluster(num_cores)
registerDoParallel(cl)

inversion_numbers <- c(1,2,3,4,5,6,7,8,9,10)
iteration_numbers <- c(1,5,10,20,50,100)
replicate_n <- 1



# Initialize an empty list to store results
all_results <- list()

# Parallelize the loop

all_results <- foreach(inversion_n = inversion_numbers, .combine = 'rbind') %:%
  foreach(iteration_n = iteration_numbers) %dopar% {

    print(paste0("Inversions:", inversion_n, ", Iterations:", iteration_n))
    tmp_result <- benchmark_sim(inversions = inversion_n, iterations = iteration_n, replicates = replicate_n)
    return(tmp_result)
  }

# Stop the parallel processing
stopCluster(cl)

all_results <- do.call(rbind, all_results)


