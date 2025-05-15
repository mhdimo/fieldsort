#ifndef FIELD_SPACE_SORT_HPP
#define FIELD_SPACE_SORT_HPP

#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <array>
#include <queue>
#include <bitset>
#include <thread>

// Fast Manhattan distance calculation - removed unnecessary template parameter
inline int manhattan_distance(const std::vector<size_t>& a, const std::vector<size_t>& b) {
    int dist = 0;
    const size_t size = a.size();
    for (size_t i = 0; i < size; ++i) {
        dist += std::abs(static_cast<int>(a[i]) - static_cast<int>(b[i]));
    }
    return dist;
}

template<typename T>
void field_space_sort_nd(std::vector<T>& data, const std::vector<size_t>& dims, size_t max_iters = 100) {
    const size_t total_elements = data.size();
    if (total_elements <= 1) return;

    // Find min and max in one pass
    auto [min_it, max_it] = std::minmax_element(data.begin(), data.end());
    T min_val = *min_it;
    T max_val = *max_it;
    
    if (min_val == max_val) return;
    
    const double range_inv = 1.0 / (max_val - min_val);
    const size_t dims_size = dims.size();

    // Pre-allocate ideals with proper size
    std::vector<std::vector<size_t>> ideals(total_elements, std::vector<size_t>(dims_size));
    
    // Pre-compute ideal positions
    for (size_t i = 0; i < total_elements; ++i) {
        double norm = static_cast<double>(data[i] - min_val) * range_inv;
        for (size_t d = 0; d < dims_size; ++d) {
            ideals[i][d] = static_cast<size_t>(norm * (dims[d] - 1));
        }
    }

    // Pre-allocate coordinate vectors
    std::vector<size_t> coord(dims_size);
    std::vector<size_t> new_coord(dims_size);
    
    // Use a priority queue instead of sorting
    using SwapPair = std::pair<int, std::pair<size_t, size_t>>;
    
    // Pre-calculate stride for faster coordinate->index conversion
    std::vector<size_t> strides(dims_size);
    strides[dims_size - 1] = 1;
    for (size_t i = dims_size - 1; i > 0; --i) {
        strides[i - 1] = strides[i] * dims[i];
    }
    
    for (size_t iter = 0; iter < max_iters; ++iter) {
        bool moved = false;
        std::vector<bool> swapped(total_elements, false);
        
        // Use a priority queue to automatically prioritize swaps with highest energy reduction
        std::priority_queue<SwapPair> swap_queue;
        
        for (size_t idx = 0; idx < total_elements; ++idx) {
            if (swapped[idx]) continue;

            // Calculate coordinates
            size_t temp = idx;
            for (size_t i = dims_size; i > 0; --i) {
                size_t d_idx = i - 1;
                coord[d_idx] = temp % dims[d_idx];
                temp /= dims[d_idx];
            }

            // Check neighbors in each dimension
            for (size_t d = 0; d < dims_size; ++d) {
                for (int delta : {-1, 1}) {
                    int coord_d = static_cast<int>(coord[d]) + delta;
                    
                    if (coord_d >= 0 && coord_d < static_cast<int>(dims[d])) {
                        // Copy current coordinate to new_coord
                        std::copy(coord.begin(), coord.end(), new_coord.begin());
                        new_coord[d] = coord_d;
                        
                        // Calculate neighbor index
                        size_t neighbor_idx = 0;
                        for (size_t i = 0; i < dims_size; ++i) {
                            neighbor_idx += new_coord[i] * strides[i];
                        }

                        if (neighbor_idx >= total_elements) continue;
                        if (swapped[neighbor_idx]) continue;

                        int current_energy = manhattan_distance(coord, ideals[idx]) +
                                           manhattan_distance(new_coord, ideals[neighbor_idx]);
                        
                        int new_energy = manhattan_distance(new_coord, ideals[idx]) +
                                      manhattan_distance(coord, ideals[neighbor_idx]);

                        int energy_reduction = current_energy - new_energy;
                        if (energy_reduction > 0) {
                            swap_queue.push({energy_reduction, {idx, neighbor_idx}});
                        }
                    }
                }
            }
        }
        
        // Process swaps in order of energy reduction (highest first)
        while (!swap_queue.empty()) {
            auto [energy_reduction, swap_pair] = swap_queue.top();
            swap_queue.pop();
            
            auto [i, j] = swap_pair;
            
            if (!swapped[i] && !swapped[j]) {
                std::swap(data[i], data[j]);
                swapped[i] = true;
                swapped[j] = true;
                moved = true;
            }
        }

        if (!moved) break;
    }
}

// 1D specialization
template<typename T>
void field_space_sort(std::vector<T>& arr, size_t max_iters = 100) {
    field_space_sort_nd(arr, {arr.size()}, max_iters);
}

// 2D specialization
template<typename T>
void field_space_sort(std::vector<std::vector<T>>& matrix, size_t max_iters = 100) {
    if (matrix.empty()) return;
    
    const size_t rows = matrix.size();
    const size_t cols = matrix[0].size();
    std::vector<T> flat;
    flat.reserve(rows * cols);
    
    for (const auto& row : matrix) {
        flat.insert(flat.end(), row.begin(), row.end());
    }
    
    field_space_sort_nd(flat, {rows, cols}, max_iters);
    
    size_t idx = 0;
    for (auto& row : matrix) {
        for (auto& val : row) {
            val = flat[idx++];
        }
    }
}

// 3D specialization
template<typename T>
void field_space_sort(std::vector<std::vector<std::vector<T>>>& cube, size_t max_iters = 100) {
    if (cube.empty() || cube[0].empty()) return;
    
    const size_t depth = cube.size();
    const size_t rows = cube[0].size();
    const size_t cols = cube[0][0].size();
    std::vector<T> flat;
    flat.reserve(depth * rows * cols);
    
    for (const auto& layer : cube) {
        for (const auto& row : layer) {
            flat.insert(flat.end(), row.begin(), row.end());
        }
    }
    
    field_space_sort_nd(flat, {depth, rows, cols}, max_iters);
    
    size_t idx = 0;
    for (auto& layer : cube) {
        for (auto& row : layer) {
            for (auto& val : row) {
                val = flat[idx++];
            }
        }
    }
}

#endif // FIELD_SPACE_SORT_HPP