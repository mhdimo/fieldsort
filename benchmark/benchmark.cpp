#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <chrono>
#include <functional>
#include "../include/field_space_sort.hpp"
#include <numeric>

template<typename Func>
long long benchmark(Func&& func) {
    auto start = std::chrono::high_resolution_clock::now();
    func();
    auto end = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
}

// 1D Benchmark
void benchmark_1d(size_t N) {
    std::vector<int> data(N);
    std::iota(data.begin(), data.end(), 0);
    std::shuffle(data.begin(), data.end(), std::mt19937{std::random_device{}()});
    
    auto original = data;
    
    long long t_std = benchmark([&] {
        std::sort(data.begin(), data.end());
    });
    
    data = original;
    long long t_field = benchmark([&] {
        field_space_sort(data);
    });

    std::cout << "1D (" << N << " elements)\n";
    std::cout << "  std::sort:  " << t_std << " ms\n";
    std::cout << "  field_sort: " << t_field << " ms\n";
    std::cout << "  Ratio: " << static_cast<float>(t_field)/t_std << "x\n\n";
}

// Fixed 2D Benchmark
void benchmark_2d(size_t rows, size_t cols) {
    const size_t total = rows * cols;
    std::vector<int> flat(total);
    std::iota(flat.begin(), flat.end(), 0);
    std::shuffle(flat.begin(), flat.end(), std::mt19937{std::random_device{}()});
    
    // Create matrix and fill properly
    std::vector<std::vector<int>> matrix(rows, std::vector<int>(cols));
    size_t idx = 0;
    for (auto& row : matrix) {
        for (auto& val : row) {
            val = flat[idx++];
        }
    }
    
    auto original = matrix;
    
    long long t_std = benchmark([&] {
        // Proper flattened sorting
        std::vector<int> sorted_flat = flat;
        std::sort(sorted_flat.begin(), sorted_flat.end());
        
        // Rebuild matrix correctly
        idx = 0;
        for (auto& row : matrix) {
            for (auto& val : row) {
                val = sorted_flat[idx++];
            }
        }
    });
    
    matrix = original;
    long long t_field = benchmark([&] {
        field_space_sort(matrix);
    });

    std::cout << "2D (" << rows << "x" << cols << ")\n";
    std::cout << "  std::sort:  " << t_std << " ms\n";
    std::cout << "  field_sort: " << t_field << " ms\n";
    std::cout << "  Ratio: " << static_cast<float>(t_field)/t_std << "x\n\n";
}

// Fixed 3D Benchmark
void benchmark_3d(size_t x, size_t y, size_t z) {
    const size_t total = x * y * z;
    std::vector<int> flat(total);
    std::iota(flat.begin(), flat.end(), 0);
    std::shuffle(flat.begin(), flat.end(), std::mt19937{std::random_device{}()});
    
    // Create cube and fill properly
    std::vector<std::vector<std::vector<int>>> cube(x, 
        std::vector<std::vector<int>>(y, std::vector<int>(z)));
    size_t idx = 0;
    for (auto& layer : cube) {
        for (auto& row : layer) {
            for (auto& val : row) {
                val = flat[idx++];
            }
        }
    }
    
    auto original = cube;
    
    long long t_std = benchmark([&] {
        // Proper flattened sorting
        std::vector<int> sorted_flat = flat;
        std::sort(sorted_flat.begin(), sorted_flat.end());
        
        // Rebuild cube correctly
        idx = 0;
        for (auto& layer : cube) {
            for (auto& row : layer) {
                for (auto& val : row) {
                    val = sorted_flat[idx++];
                }
            }
        }
    });
    
    cube = original;
    long long t_field = benchmark([&] {
        field_space_sort(cube);
    });

    std::cout << "3D (" << x << "x" << y << "x" << z << ")\n";
    std::cout << "  std::sort:  " << t_std << " ms\n";
    std::cout << "  field_sort: " << t_field << " ms\n";
    std::cout << "  Ratio: " << static_cast<float>(t_field)/t_std << "x\n\n";
}

int main() {
    // 1D Tests
    benchmark_1d(10'000);
    benchmark_1d(50'000);

    // 2D Tests
    benchmark_2d(100, 100);  // 10,000 elements
    benchmark_2d(50, 200);   // 10,000 elements

    // 3D Tests
    benchmark_3d(20, 20, 20);  // 8,000 elements
    benchmark_3d(30, 10, 10);  // 3,000 elements

    return 0;
}