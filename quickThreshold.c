#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>

#define DEFAULT_THRESHOLD 25 // Default threshold for switching to insertion sort

// Insertion sort function for small sub-arrays
void insertion_sort(int *arr, int left, int right) {
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

// Function to partition the array
int partition(int *arr, int low, int high) {
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

// Quicksort function with insertion sort for small sub-arrays
void quicksort(int *arr, int low, int high, int threshold) {
    if (low < high) {
        if (high - low < threshold) {
            insertion_sort(arr, low, high);
        } else {
            int pi = partition(arr, low, high);
            quicksort(arr, low, pi - 1, threshold);
            quicksort(arr, pi + 1, high, threshold);
        }
    }
}

// Min-heap priority queue node structure
typedef struct {
    int value;
    int array_idx;  // Index of the array from which the value came
} MinHeapNode;

// Swap two heap nodes
void swap(MinHeapNode *x, MinHeapNode *y) {
    MinHeapNode temp = *x;
    *x = *y;
    *y = temp;
}

// Heapify the heap
void heapify(MinHeapNode *heap, int size, int i) {
    int smallest = i;
    int left = 2 * i + 1;
    int right = 2 * i + 2;

    if (left < size && heap[left].value < heap[smallest].value) {
        smallest = left;
    }

    if (right < size && heap[right].value < heap[smallest].value) {
        smallest = right;
    }

    if (smallest != i) {
        swap(&heap[i], &heap[smallest]);
        heapify(heap, size, smallest);
    }
}

// K-Way Merge function using a min-heap priority queue
void k_way_merge(int *arr, int size, int num_procs, int chunk_size) {
    int *indices = (int *)malloc(num_procs * sizeof(int));  // Track indices of each sub-array
    int *output = (int *)malloc(size * sizeof(int));  // Array to store the final sorted result
    MinHeapNode *heap = (MinHeapNode *)malloc(num_procs * sizeof(MinHeapNode));

    // Initialize the heap with the first element of each sub-array
    for (int i = 0; i < num_procs; i++) {
        indices[i] = i * chunk_size;  // Set starting index of each sub-array
        if (indices[i] < size) {
            heap[i].value = arr[indices[i]];
            heap[i].array_idx = i;
        } else {
            heap[i].value = INT_MAX;  // Mark empty sub-array
        }
    }

    // Build the initial min-heap
    for (int i = (num_procs - 1) / 2; i >= 0; i--) {
        heapify(heap, num_procs, i);
    }

    // Perform the merge
    for (int i = 0; i < size; i++) {
        // Extract the smallest element from the heap
        MinHeapNode root = heap[0];
        output[i] = root.value;

        // Get the next element from the corresponding sub-array
        indices[root.array_idx]++;
        if (indices[root.array_idx] < (root.array_idx + 1) * chunk_size && indices[root.array_idx] < size) {
            heap[0].value = arr[indices[root.array_idx]];
        } else {
            heap[0].value = INT_MAX;  // Mark the sub-array as exhausted
        }

        // Heapify the root element
        heapify(heap, num_procs, 0);
    }

    // Copy the merged result back to the original array
    for (int i = 0; i < size; i++) {
        arr[i] = output[i];
    }

    free(indices);
    free(output);
    free(heap);
}

int main(int argc, char *argv[]) {
    int rank, size, N, padded_N;
    int *arr = NULL;
    double start_time, end_time, total_time, process_time;
    
    // Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Set a fixed seed for reproducibility
    srand(42); // Fixed seed

    // Get the number of elements to sort from the user
    if (rank == 0) {
        printf("Enter the number of elements to sort: ");
        scanf("%d", &N);
    }

    // Broadcast N to all processes
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Calculate chunk size and padded size
    int chunk_size = (N + size - 1) / size;  // Ensure enough space
    padded_N = chunk_size * size;

    // Allocate memory for the array on the root process
    if (rank == 0) {
        arr = (int *)malloc(padded_N * sizeof(int));
        // Generate random input on the root process
        for (int i = 0; i < N; i++) {
            arr[i] = rand() % 1000000;  // Random values between 0 and 999999
        }
        // Padding the remaining elements (if any)
        for (int i = N; i < padded_N; i++) {
            arr[i] = INT_MAX;  // Fill extra space with large values
        }
    }

    // Start timing the overall program
    start_time = MPI_Wtime();

    // Allocate space for local array
    int *local_arr = (int *)malloc(chunk_size * sizeof(int));

    // Scatter the array to all processes
    MPI_Scatter(arr, chunk_size, MPI_INT, local_arr, chunk_size, MPI_INT, 0, MPI_COMM_WORLD);

    // Start timing the process-specific time
    process_time = MPI_Wtime();

    // Check for a threshold value from command line argument or use default
    int threshold = DEFAULT_THRESHOLD;
    if (argc > 1) {
        threshold = atoi(argv[1]);
    }

    // Apply parallel quicksort with threshold check
    quicksort(local_arr, 0, chunk_size - 1, threshold);

    // End process-specific time
    process_time = MPI_Wtime() - process_time;

    // Gather the sorted arrays from all processes
    MPI_Gather(local_arr, chunk_size, MPI_INT, arr, chunk_size, MPI_INT, 0, MPI_COMM_WORLD);

    // Merge sorted sub-arrays using K-Way Merge in the root process
    if (rank == 0) {
        k_way_merge(arr, N, size, chunk_size);
        printf("Final sorted array: ");
        for (int i = 0; i < N; i++) {
            printf("%d ", arr[i]);
        }
        printf("\n");
    }

    // End timing and calculate the total execution time
    end_time = MPI_Wtime();
    total_time = end_time - start_time;

    // Print execution time in the root process
    if (rank == 0) {
        printf("Program execution time: %f seconds\n", total_time);
    }

    // Print process-specific time
    printf("Process %d execution time: %f seconds\n", rank, process_time);

    // Free allocated memory
    free(local_arr);

    if (rank == 0) {
        free(arr);
    }

    // Finalize MPI
    MPI_Finalize();
    return 0;
}
