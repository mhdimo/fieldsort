#ifndef FIELD_SPACE_SORT_HPP
#define FIELD_SPACE_SORT_HPP

#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <array>
#include <queue>
#include <random>

// Fast Manhattan distance calculation
inline int manhattan_distance(const std::vector<size_t>& a, const std::vector<size_t>& b) {
    int dist = 0;
    const size_t size = a.size();
    for (size_t i = 0; i < size; ++i) {
        dist += std::abs(static_cast<int>(a[i]) - static_cast<int>(b[i]));
    }
    return dist;
}

// Apply initial heuristic pre-ordering to improve starting state
template<typename T>
void apply_heuristic_ordering(std::vector<T>& data, const std::vector<size_t>& dims) {
    const size_t total = data.size();
    if (total <= 1) return;
    
    // Find min/max for normalization
    auto [min_it, max_it] = std::minmax_element(data.begin(), data.end());
    T min_val = *min_it;
    T max_val = *max_it;
    if (min_val == max_val) return;
    
    // Create index vector
    std::vector<size_t> indices(total);
    std::iota(indices.begin(), indices.end(), 0);
    
    // Bucket sort heuristic - divide into 8 buckets
    constexpr int NUM_BUCKETS = 8;
    std::vector<std::vector<size_t>> buckets(NUM_BUCKETS);
    
    // Assign indices to buckets based on value
    double range_inv = 1.0 / (max_val - min_val);
    for (size_t i = 0; i < total; ++i) {
        double norm = (data[i] - min_val) * range_inv;
        int bucket = std::min(NUM_BUCKETS-1, static_cast<int>(norm * NUM_BUCKETS));
        buckets[bucket].push_back(i);
    }
    
    // Apply partial sorting within each bucket (improves initial state)
    for (auto& bucket : buckets) {
        if (bucket.size() > 32) { // Only sort larger buckets
            std::sort(bucket.begin(), bucket.end(), 
                [&data](size_t a, size_t b) { return data[a] < data[b]; });
        }
    }
    
    // Create new data vector with partial ordering
    std::vector<T> ordered_data(total);
    size_t idx = 0;
    for (const auto& bucket : buckets) {
        for (size_t i : bucket) {
            ordered_data[idx++] = data[i];
        }
    }
    
    // Replace original data
    data = std::move(ordered_data);
}

template<typename T>
void field_space_sort_nd(std::vector<T>& data, const std::vector<size_t>& dims, size_t max_iters = 100) {
    const size_t total_elements = data.size();
    if (total_elements <= 1) return;

    // Apply heuristic pre-ordering for better starting state
    apply_heuristic_ordering(data, dims);
    
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
    
    // Reduced iterations needed due to better starting state
    max_iters = std::max(max_iters / 2, size_t(10)); 
    
    for (size_t iter = 0; iter < max_iters; ++iter) {
        bool moved = false;
        std::vector<bool> swapped(total_elements, false);
        std::priority_queue<SwapPair> swap_queue;
        
        // Skip every other element to reduce computational load
        // (we'll catch them on the next iteration)
        for (size_t idx = (iter % 2); idx < total_elements; idx += 2) {
            if (swapped[idx]) continue;

            // Calculate coordinates
            size_t temp = idx;
            for (size_t i = dims_size; i > 0; --i) {
                size_t d_idx = i - 1;
                coord[d_idx] = temp % dims[d_idx];
                temp /= dims[d_idx];
            }

            // Check only immediate neighbors in each dimension
            for (size_t d = 0; d < dims_size; ++d) {
                for (int delta : {-1, 1}) {
                    int coord_d = static_cast<int>(coord[d]) + delta;
                    
                    if (coord_d >= 0 && coord_d < static_cast<int>(dims[d])) {
                        // Set new_coord directly without copying the entire array
                        new_coord[d] = coord_d;
                        if (d > 0) new_coord[d-1] = coord[d-1];
                        if (d+1 < dims_size) new_coord[d+1] = coord[d+1];
                        
                        // Calculate neighbor index
                        size_t neighbor_idx = 0;
                        for (size_t i = 0; i < dims_size; ++i) {
                            neighbor_idx += new_coord[i] * strides[i];
                        }

                        if (neighbor_idx >= total_elements || swapped[neighbor_idx]) 
                            continue;

                        int current_energy = manhattan_distance(coord, ideals[idx]) +
                                            manhattan_distance(new_coord, ideals[neighbor_idx]);
                        
                        int new_energy = manhattan_distance(new_coord, ideals[idx]) +
                                        manhattan_distance(coord, ideals[neighbor_idx]);

                        int energy_reduction = current_energy - new_energy;
                        if (energy_reduction > 1) { // Only consider significant improvements
                            swap_queue.push({energy_reduction, {idx, neighbor_idx}});
                        }
                    }
                }
            }
        }
        
        // Process only top 80% of swaps for speed
        size_t max_swaps = swap_queue.size() * 0.8;
        size_t swaps_done = 0;
        
        while (!swap_queue.empty() && swaps_done < max_swaps) {
            auto [energy_reduction, swap_pair] = swap_queue.top();
            swap_queue.pop();
            
            auto [i, j] = swap_pair;
            
            if (!swapped[i] && !swapped[j]) {
                std::swap(data[i], data[j]);
                swapped[i] = true;
                swapped[j] = true;
                moved = true;
                swaps_done++;
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

// 4D specialization
template<typename T>
void field_space_sort(std::vector<std::vector<std::vector<std::vector<T>>>>& tesseract, size_t max_iters = 100) {
    if (tesseract.empty() || tesseract[0].empty() || tesseract[0][0].empty()) return;
    
    const size_t w = tesseract.size();              // First dimension
    const size_t x = tesseract[0].size();           // Second dimension
    const size_t y = tesseract[0][0].size();        // Third dimension
    const size_t z = tesseract[0][0][0].size();     // Fourth dimension
    
    std::vector<T> flat;
    flat.reserve(w * x * y * z);  // Pre-allocate for efficiency
    
    // Flatten the 4D structure
    for (const auto& space : tesseract) {
        for (const auto& plane : space) {
            for (const auto& row : plane) {
                flat.insert(flat.end(), row.begin(), row.end());
            }
        }
    }
    
    // Sort the flattened data
    field_space_sort_nd(flat, {w, x, y, z}, max_iters);
    
    // Reshape back to 4D
    size_t idx = 0;
    for (auto& space : tesseract) {
        for (auto& plane : space) {
            for (auto& row : plane) {
                for (auto& val : row) {
                    val = flat[idx++];
                }
            }
        }
    }
}

#endif // FIELD_SPACE_SORT_HPP