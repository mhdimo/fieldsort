# FieldSort: A C++ Library for N-Dimensional Data Sorting

FieldSort is a C++ header-only library that provides a custom sorting algorithm, `field_space_sort`, designed for sorting elements in N-dimensional data structures. The algorithm aims to arrange data based on proximity and a concept of "energy reduction" within the N-dimensional space.

## Features

*   **N-Dimensional Sorting**: Supports sorting for:
    *   1D `std::vector<T>`
    *   2D `std::vector<std::vector<T>>` (matrices)
    *   3D `std::vector<std::vector<std::vector<T>>>` (cubes)
    *   4D `std::vector<std::vector<std::vector<std::vector<T>>>>` (tesseracts)
    *   A generic `field_space_sort_nd` for flat arrays with dimension information.
*   **Heuristic Pre-ordering**: Applies an initial heuristic ordering to potentially improve the starting state for the main sorting algorithm.
*   **Multiple Algorithm Versions**:
    *   `field_space_sort.hpp`: The original implementation.
    *   `field_space_sortv2.hpp`: An alternative version with different heuristics and swap processing.
    *   `field_space_sort_parallel.hpp`: A version incorporating parallelism for potentially faster execution on multi-core processors.
*   **Templated**: Works with various data types `T` that support comparison and basic arithmetic.
*   **Benchmarking Suite**: Includes benchmark programs to compare the performance of `field_space_sort` against `std::sort`.
*   **Test Cases**: Provides basic test cases to verify the functionality of the sorting algorithm.

## Getting Started

### Prerequisites

*   A C++ compiler that supports C++17 (due to usage of `std::shuffle`, `std::iota`, structured bindings, etc.).
*   Standard C++ libraries (vector, algorithm, cmath, numeric, chrono, etc.).
*   (Optional, for parallel version) Threading support (`<thread>`, `<mutex>`, `<future>`).

### Building

FieldSort is a header-only library, so there's no separate compilation step for the library itself. You just need to include the desired header file in your project.

To compile the provided test cases or benchmark programs, you can use a C++ compiler like g++.

**Example: Compiling Test Cases (using `field_space_sort.hpp`)**

```bash
g++ -std=c++17 -I./include test/test_cases.cpp -o test_cases
```

**Example: Compiling a Benchmark (e.g., `benchmark.cpp`)**

```bash
g++ -std=c++17 -I./include benchmark/benchmark.cpp -o benchmark -O3 # -O3 for performance
```

Replace `benchmark.cpp` with `benchmarkv2.cpp` or `benchmarkparallel.cpp` and the corresponding header in the include path if needed (though the benchmarks include the correct headers relative to their location). For `benchmarkparallel.cpp`, you might need to link pthreads:

```bash
g++ -std=c++17 -I./include benchmark/benchmarkparallel.cpp -o benchmarkparallel -O3 -pthread
```

### Running Tests

After compiling the test cases:

```bash
./test_cases
```

You should see "All tests passed!" if everything is correct.

### Running Benchmarks

After compiling a benchmark program:

```bash
./benchmark
# or ./benchmarkv2
# or ./benchmarkparallel
```

The output will show the time taken by `std::sort` versus `field_space_sort` for various data sizes and dimensions, along with a ratio.

## Usage

To use `field_space_sort` in your project, include the relevant header file.

### Include Header

Choose one of the following based on the version you want to use:

```cpp
#include "field_space_sort.hpp" // For the original version
// or
#include "field_space_sortv2.hpp" // For the V2 version
// or
#include "field_space_sort_parallel.hpp" // For the parallel version
```

### Basic Example (1D Vector)

```cpp
#include <iostream>
#include <vector>
#include "field_space_sort.hpp" // Or your chosen version

int main() {
    std::vector<int> arr = {9, 5, 3, 7, 1, 8, 2, 6, 4, 0};

    std::cout << "Original array: ";
    for (int x : arr) {
        std::cout << x << " ";
    }
    std::cout << std::endl;

    field_space_sort(arr); // Default max_iters = 100

    std::cout << "Sorted array:   ";
    for (int x : arr) {
        std::cout << x << " ";
    }
    std::cout << std::endl;

    return 0;
}
```

### Example (2D Matrix)

```cpp
#include <iostream>
#include <vector>
#include "field_space_sort.hpp" // Or your chosen version

int main() {
    std::vector<std::vector<int>> matrix = {
        {9, 5, 3},
        {7, 1, 8},
        {2, 6, 4}
    };

    std::cout << "Original matrix:" << std::endl;
    for (const auto& row : matrix) {
        for (int x : row) {
            std::cout << x << " ";
        }
        std::cout << std::endl;
    }

    field_space_sort(matrix); // Default max_iters = 100

    std::cout << "\nSorted matrix:" << std::endl;
    for (const auto& row : matrix) {
        for (int x : row) {
            std::cout << x << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}
```

### Using `field_space_sort_nd` (N-Dimensional Flat Array)

The `_nd` versions of the sort function (e.g., `field_space_sort_nd` in [`include/field_space_sort.hpp`](include/field_space_sort.hpp), [`include/field_space_sortv2.hpp`](include/field_space_sortv2.hpp), and [`include/field_space_sort_parallel.hpp`](include/field_space_sort_parallel.hpp)) take a flat `std::vector<T>` and a `std::vector<size_t>` representing the dimensions.

```cpp
#include <iostream>
#include <vector>
#include <numeric> // For std::iota
#include <algorithm> // For std::shuffle
#include <random> // For std::mt19937
#include "field_space_sort.hpp" // Or your chosen version

int main() {
    std::vector<size_t> dims = {2, 3, 2}; // Example 2x3x2 structure
    size_t total_elements = 1;
    for (size_t d : dims) {
        total_elements *= d;
    }

    std::vector<int> data(total_elements);
    std::iota(data.begin(), data.end(), 0); // Fill with 0, 1, 2, ...
    std::shuffle(data.begin(), data.end(), std::mt19937{std::random_device{}()}); // Shuffle

    std::cout << "Original flat data: ";
    for (int x : data) {
        std::cout << x << " ";
    }
    std::cout << std::endl;

    field_space_sort_nd(data, dims);

    std::cout << "Sorted flat data:   ";
    for (int x : data) {
        std::cout << x << " ";
    }
    std::cout << std::endl;

    return 0;
}
```

The `max_iters` parameter can be passed to all `field_space_sort` and `field_space_sort_nd` functions to control the maximum number of iterations the algorithm performs.

## Algorithm Versions

*   **`include/field_space_sort.hpp`**: This is the original implementation of the FieldSort algorithm. It uses Manhattan distance to calculate "energy" and iteratively swaps elements to reduce this energy.
*   **`include/field_space_sortv2.hpp`**: This version introduces some modifications, such as:
    *   A different heuristic pre-ordering step ([`apply_heuristic_ordering`](include/field_space_sortv2.hpp#L45)).
    *   Processing only a percentage of potential swaps ([`field_space_sort_nd`](include/field_space_sortv2.hpp#L172)).
    *   Considering only swaps with an energy reduction greater than 1 ([`field_space_sort_nd`](include/field_space_sortv2.hpp#L148)).
*   **`include/field_space_sort_parallel.hpp`**: This version aims to leverage multi-core processors by parallelizing parts of the sorting algorithm, particularly the heuristic pre-ordering and potentially other loops within the main sort. It includes functions like [`get_optimal_thread_count`](include/field_space_sort_parallel.hpp#L21).

## Directory Structure

```
fieldsort/
├── include/
│   ├── field_space_sort.hpp         # Original version
│   ├── field_space_sort_parallel.hpp # Parallel version
│   └── field_space_sortv2.hpp       # V2 version
├── benchmark/
│   ├── benchmark.cpp                # Benchmarks for field_space_sort.hpp
│   ├── benchmarkparallel.cpp        # Benchmarks for field_space_sort_parallel.hpp
│   └── benchmarkv2.cpp              # Benchmarks for field_space_sortv2.hpp
├── test/
│   └── test_cases.cpp               # Test cases (uses field_space_sort.hpp by default)
└── README.md                        # This file
```

## Contributing

Contributions are welcome! Please feel free to submit pull requests or open issues for bugs, feature requests, or improvements.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
