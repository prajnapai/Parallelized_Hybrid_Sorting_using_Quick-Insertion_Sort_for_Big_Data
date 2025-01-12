#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define THRESHOLD 25  // Threshold for switching to Insertion Sort

// Function to perform Insertion Sort
void insertionSort(int arr[], int left, int right) {
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

// Function to partition the array for Quick Sort
int partition(int arr[], int low, int high) {
    int pivot = arr[high];
    int i = low - 1;
    for (int j = low; j < high; j++) {
        if (arr[j] < pivot) 
        {
            i++;
            int temp = arr[i];
            arr[i] = arr[j];
            arr[j] = temp;
        }
    }
    int temp = arr[i + 1];
    arr[i + 1] = arr[high];
    arr[high] = temp;
    return i + 1;
}

// Hybrid Quick Sort with Insertion Sort
void quickSort(int arr[], int low, int high) {
    if (low < high) {
        // If the sub-array size is below the threshold, use Insertion Sort
        if (high - low  < THRESHOLD) {
            insertionSort(arr, low, high);
        } else {
            // Otherwise, use Quick Sort
            int pivotIndex = partition(arr, low, high);
            quickSort(arr, low, pivotIndex - 1);
            quickSort(arr, pivotIndex + 1, high);
        }
    }
}

int main() {
    int N;
    printf("Enter the number of elements to sort:\n");
    scanf("%d", &N);

    // Dynamically allocate array and fill it with random numbers
    int *arr = (int *)malloc(N * sizeof(int));
    if (arr == NULL) {
        printf("Memory allocation failed!\n");
        return -1;
    }

    srand(time(0));
    for (int i = 0; i < N; i++) {
            arr[i] = ((long long)rand() * rand()) % 5000000;
        // arr[i] = rand() % 1000000;  // Random numbers between 0 and 999999
    }

    printf("Unsorted array: \n");
    for (int i = 0; i < N; i++) {
        printf("%d ", arr[i]);
    }
    printf("\n");

    // Measure the time taken for sorting
    clock_t start, end;

    start = clock();
    // Perform Quick Sort with Insertion Sort for small sub-arrays
    quickSort(arr, 0, N - 1);
    end = clock();

    double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;

    // Print the sorted array
    printf("Sorted array: \n");
    for (int i = 0; i < N; i++) {
        printf("%d ", arr[i]);
    }
    printf("\n");

    printf("Time taken for sorting: %f seconds\n", time_taken);

    // Free dynamically allocated memory
    free(arr);

    return 0;
}
