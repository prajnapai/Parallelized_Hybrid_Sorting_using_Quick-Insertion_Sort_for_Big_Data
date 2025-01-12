Project Overview

This project implements a Parallelized Hybrid Sorting Algorithm combining Quick Sort and Insertion Sort for efficient handling of large datasets. The hybrid approach leverages the divide-and-conquer technique of Quick Sort with the simplicity and efficiency of Insertion Sort for smaller partitions. Parallelization using MPI and CUDA further optimizes the sorting process for big data applications.

Methodology
1. Data: Randomly generated numbers were used for various inputs ranging from 100 elements to 4000000 elements for performance testing.
2. Hybrid Sorting:
   * Quick Sort divides the dataset into partitions.
   * Insertion Sort handles smaller partitions for improved performance when a certain threshold is reached.
3. Parallel Implementation:
   * MPI: Divides the dataset among multiple processors.
   * CUDA: Parallelizes the sorting on GPUs.
  

Results and Performance
Datasets Tested: 1000, 50000, 100000, 1000000, 2000000, 4000000 elements.
Processors: 4, 8, 16, 24.
Performance Metrics: Execution time, speedup, efficiency.
Observations: MPI showed better performance for smaller datasets, while CUDA excelled for massive data sizes.
