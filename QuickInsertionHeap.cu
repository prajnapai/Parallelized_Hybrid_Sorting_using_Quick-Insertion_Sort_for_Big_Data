#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>

#define THRESHOLD 25  // Threshold for switching to Insertion Sort
#define MAX_STACK_SIZE 4096  // Maximum stack size for iterative quicksort

// Insertion Sort (Device Function)
__device__ void insertionSort(int *arr, int left, int right) {
    for (int i = left + 1; i <= right; i++) {
        int key = arr[i];
        int j = i - 1;
        while (j >= left && arr[j] > key) {
            arr[j + 1] = arr[j];
            j--;
        }
        arr[j + 1] = key;
    }
}

// Partition Function for Quick Sort (Device Function)
__device__ int partition(int *arr, int low, int high) {
    int pivot = arr[high];
    int i = (low - 1);
    for (int j = low; j < high; j++) {
        if (arr[j] < pivot) {
            i++;
            int temp = arr[i];
            arr[i] = arr[j];
            arr[j] = temp;
        }
    }
    int temp = arr[i + 1];
    arr[i + 1] = arr[high];
    arr[high] = temp;
    return (i + 1);
}

// Iterative Quick Sort (Device Function)
__device__ void iterativeQuickSort(int *arr, int low, int high) {
    if (low >= high) return; // Avoid unnecessary sorting

    int stack[MAX_STACK_SIZE];
    int top = -1;

    stack[++top] = low;
    stack[++top] = high;

    while (top >= 0) {
        high = stack[top--];
        low = stack[top--];

        // If the size of the current partition is small, use insertion sort
        if (high - low < THRESHOLD) {
            insertionSort(arr, low, high);
            continue;  // Skip to the next iteration
        }

        int pi = partition(arr, low, high);

        // Push left side to stack
        if (pi - 1 > low) {
            stack[++top] = low;
            stack[++top] = pi - 1;
        }
        // Push right side to stack
        if (pi + 1 < high) {
            stack[++top] = pi + 1;
            stack[++top] = high;
        }
    }
}

// Kernel for Parallel Quick Sort
__global__ void parallelQuickSort(int *arr, int size, int segmentSize) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    // Each thread works on a segment of the array
    int start = idx * segmentSize;
    int end = min(start + segmentSize - 1, size - 1);

    // Ensure valid indices before sorting
    if (start <= end) {
        iterativeQuickSort(arr, start, end);
    }
}

// Kernel function to perform K-Way Merge using Min-Heap
__global__ void kWayMerge(int *arr, int *temp, int size, int segmentSize) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    // Each thread merges its own segments based on K-Way merging
    int start = idx * segmentSize * 2;
    int end1 = start + segmentSize - 1;
    int start2 = end1 + 1;
    int end2 = start2 + segmentSize - 1;

    if (start < size) {
        if (end1 >= size) end1 = size - 1;
        
        if (start2 < size) {
            if (end2 >= size) end2 = size - 1;

            int i = start, j = start2, k = start;
            while (i <= end1 && j <= end2) {
                temp[k++] = (arr[i] <= arr[j]) ? arr[i++] : arr[j++];
            }
            while (i <= end1) {
                temp[k++] = arr[i++];
            }
            while (j <= end2) {
                temp[k++] = arr[j++];
            }
        }
    }
}

void writeArrayToFile(int *arr, int size) {
    FILE *file = fopen("sorted.txt", "w");
    if (file == NULL) {
        printf("Error opening file!\n");
        return;
    }

    for (int i = 0; i < size; i++) {
        fprintf(file, "%d\n", arr[i]);
    }

    fclose(file);
    printf("Sorted output written to sorted.txt\n");
}

int main() {
    int SIZE;
    printf("Enter the number of elements: ");
    scanf("%d", &SIZE);

    int *h_arr = (int *)malloc(SIZE * sizeof(int));
    int *d_arr, *d_temp;

    // Initialize the array with random values
    srand(time(0));  // Seed for random number generation
    for (int i = 0; i < SIZE; i++) {
        h_arr[i] = rand() % 5000000;
    }

    // Print the unsorted array
    printf("Unsorted array:\n");
    for (int i = 0; i < SIZE; i++) {
        printf("%d ", h_arr[i]);
    }
    printf("\n");

    // Allocate device memory with error checking
    cudaError_t err;
    
    err = cudaMalloc((void **)&d_arr, SIZE * sizeof(int));
    if (err != cudaSuccess) {
        fprintf(stderr, "Error allocating device memory for d_arr: %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    err = cudaMalloc((void **)&d_temp, SIZE * sizeof(int));
    if (err != cudaSuccess) {
        fprintf(stderr, "Error allocating device memory for d_temp: %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    // Copy host array to device
    cudaMemcpy(d_arr, h_arr, SIZE * sizeof(int), cudaMemcpyHostToDevice);

    // Timing variables
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

   // Kernel configuration
   int threadsPerBlock = 256;
   int blocksPerGrid = (SIZE + threadsPerBlock - 1) / threadsPerBlock;

   // Determine segment size for each thread to work on
   int segmentSize = (SIZE + blocksPerGrid * threadsPerBlock - 1) / (blocksPerGrid * threadsPerBlock);

   cudaEventRecord(start, 0);
   // Parallel Quick Sort
   parallelQuickSort<<<blocksPerGrid, threadsPerBlock>>>(d_arr, SIZE, segmentSize);
   cudaDeviceSynchronize();  // Ensure sorting is complete

   // Stop timing for Quick Sort
   cudaEventRecord(stop, 0);
   cudaEventSynchronize(stop);

   // Calculate the elapsed time for Quick Sort only
   float milliseconds = 0;
   cudaEventElapsedTime(&milliseconds, start, stop);
   float seconds = milliseconds / 1000;

   // Merging sorted segments with doubling segment size using K-Way Merge
   int currentSegmentSize = segmentSize;
   while (currentSegmentSize < SIZE) {
       int mergeBlocks = (SIZE + 2 * currentSegmentSize - 1) / (2 * currentSegmentSize);
       kWayMerge<<<mergeBlocks, threadsPerBlock>>>(d_arr, d_temp, SIZE, currentSegmentSize);
       cudaMemcpy(d_arr, d_temp, SIZE * sizeof(int), cudaMemcpyDeviceToDevice);
       cudaDeviceSynchronize(); // Ensure all threads have finished before copying back
       currentSegmentSize *= 2;  // Double the segment size for the next merging step
   }

   // Copy sorted array back to host
   cudaMemcpy(h_arr, d_arr, SIZE * sizeof(int), cudaMemcpyDeviceToHost);

   // Print sorted array
   printf("Sorted array:\n");
   for (int i = 0; i < SIZE; i++) {
       printf("%d ", h_arr[i]);
   }
   printf("\n");

   printf("Total Kernel execution time: %f seconds\n", seconds);

   writeArrayToFile(h_arr, SIZE);

   // Free memory
   free(h_arr);
   cudaFree(d_arr);
   cudaFree(d_temp);
   cudaEventDestroy(start);
   cudaEventDestroy(stop);

   return 0;
}
