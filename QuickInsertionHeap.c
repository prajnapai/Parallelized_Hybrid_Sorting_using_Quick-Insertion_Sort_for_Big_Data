#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>

#define THRESHOLD 25 // Threshold for switching to insertion sort

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
void quicksort(int *arr, int low, int high) {
    if (low < high) {
        if (high - low < THRESHOLD) {
            insertion_sort(arr, low, high);
        } else {
            int pi = partition(arr, low, high);
            quicksort(arr, low, pi - 1);
            quicksort(arr, pi + 1, high);
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
void k_way_merge(int *arr, int N, int num_procs, int *sendcounts, int *displs) {
    int *indices = (int *)malloc(num_procs * sizeof(int));  // Track indices of each sub-array
    int *output = (int *)malloc(N * sizeof(int));  // Array to store the final sorted result
    MinHeapNode *heap = (MinHeapNode *)malloc(num_procs * sizeof(MinHeapNode));

    // Initialize the heap with the first element of each sub-array
    for (int i = 0; i < num_procs; i++) {
        indices[i] = displs[i];  // Set starting index of each sub-array
        if (sendcounts[i] > 0) {
            heap[i].value = arr[indices[i]];
            heap[i].array_idx = i;
            indices[i]++;  // Move to the next element in the sub-array
        } else {
            heap[i].value = INT_MAX;  // Mark empty sub-array
            heap[i].array_idx = i;
        }
    }

    // Build the initial min-heap
    for (int i = (num_procs - 1) / 2; i >= 0; i--) {
        heapify(heap, num_procs, i);
    }

    // Perform the merge
    for (int i = 0; i < N; i++) {
        // Extract the smallest element from the heap
        MinHeapNode root = heap[0];
        output[i] = root.value;

        // Get the next element from the corresponding sub-array
        if (indices[root.array_idx] < displs[root.array_idx] + sendcounts[root.array_idx]) {
            heap[0].value = arr[indices[root.array_idx]];
            indices[root.array_idx]++;
        } else {
            heap[0].value = INT_MAX;  // Mark the sub-array as exhausted
        }

        // Heapify the root element to maintain the min-heap property
        heapify(heap, num_procs, 0);
    }

    // Copy the merged result back to the original array
    for (int i = 0; i < N; i++) {
        arr[i] = output[i];
    }

    free(indices);
    free(output);
    free(heap);
}

int main(int argc, char *argv[]) {
    int rank, size, N;
    int *arr = NULL;
    double start_time, end_time, total_time;
    int *sendcounts = NULL;
    int *displs = NULL;

    // Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Seed the random number generator
    srand(time(NULL) + rank);

    // Get the number of elements to sort from the user
    if (rank == 0) {
        printf("Enter the number of elements to sort:\n");
        fflush(stdout); // Ensure prompt is printed before scanf
        scanf("%d", &N);
    }

    // Broadcast N to all processes
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Calculate sendcounts and displacements on the root process
    if (rank == 0) {
        sendcounts = (int *)malloc(size * sizeof(int));
        displs = (int *)malloc(size * sizeof(int));

        int base = N / size;
        int remainder = N % size;

        for (int i = 0; i < size; i++) {
            sendcounts[i] = base + (i < remainder ? 1 : 0);
        }

        displs[0] = 0;
        for (int i = 1; i < size; i++) {
            displs[i] = displs[i - 1] + sendcounts[i - 1];
        }

        // Allocate memory for the array on the root process
        arr = (int *)malloc(N * sizeof(int));
        // Generate random input on the root process
        for (int i = 0; i < N; i++) {
            arr[i] = ((long long)rand() * rand()) % 5000000;
            // arr[i] = rand() % 1000000;  // Random values between 0 and 999999
        }
    }

    // Start timing the overall program

    // Determine the number of elements for this process
    int local_N;
    if (rank == 0) {
        local_N = sendcounts[0];
    }
    // Scatter the sendcounts to all processes to know their local_N
    MPI_Scatter(sendcounts, 1, MPI_INT, &local_N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Allocate space for local array
    int *local_arr = (int *)malloc(local_N * sizeof(int));

    // Scatter the array to all processes using Scatterv
    MPI_Scatterv(arr, sendcounts, displs, MPI_INT, local_arr, local_N, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD); // Ensure all processes start timing together
    start_time = MPI_Wtime();
    // Apply parallel quicksort with threshold check
    quicksort(local_arr, 0, local_N - 1);
    end_time = MPI_Wtime();
    // Gather the sorted arrays from all processes using Gatherv
    MPI_Gatherv(local_arr, local_N, MPI_INT, arr, sendcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);

    // Merge sorted sub-arrays using K-Way Merge in the root process
    if (rank == 0) {
        k_way_merge(arr, N, size, sendcounts, displs);
        printf("Final sorted array: ");
        for (int i = 0; i < N; i++) {
            printf("%d ", arr[i]);
        }
        printf("\n");
    }

    MPI_Barrier(MPI_COMM_WORLD); // Ensure all processes start timing together
    // End timing and calculate the total execution time
    // end_time = MPI_Wtime();
    total_time = end_time - start_time;

    // Print execution time in the root process
    if (rank == 0) {
        printf("Program execution time: %f seconds\n", total_time);
    }

    // Free allocated memory
    free(local_arr);
    if (rank == 0) {
        free(arr);
        free(sendcounts);
        free(displs);
    }

    // Finalize MPI
    MPI_Finalize();
    return 0;
}
