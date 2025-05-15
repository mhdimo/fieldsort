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
#include <thread>
#include <mutex>
#include <atomic>
#include <future>

// Fast Manhattan distance calculation
inline int manhattan_distance(const std::vector<size_t>& a, const std::vector<size_t>& b) {
    int dist = 0;
    const size_t size = a.size();
    for (size_t i = 0; i < size; ++i) {
        dist += std::abs(static_cast<int>(a[i]) - static_cast<int>(b[i]));
    }
    return dist;
}

// Get optimal number of threads based on hardware and problem size
inline size_t get_optimal_thread_count(size_t problem_size) {
    const size_t hw_threads = std::thread::hardware_concurrency();
    const size_t min_items_per_thread = 1000; // Avoid excessive threading overhead for small problems
    
    return std::max(size_t(1), std::min(hw_threads, 
           std::max(size_t(1), problem_size / min_items_per_thread)));
}

// Apply initial heuristic pre-ordering to improve starting state - parallel version
template<typename T>
void apply_heuristic_ordering(std::vector<T>& data, const std::vector<size_t>& dims) {
    const size_t total = data.size();
    if (total <= 1) return;
    
    // Find min/max for normalization
    auto [min_it, max_it] = std::minmax_element(data.begin(), data.end());
    T min_val = *min_it;
    T max_val = *max_it;
    if (min_val == max_val) return;
    
    // Bucket sort heuristic - divide into buckets
    constexpr int NUM_BUCKETS = 16; // Increased from 8 for better parallelism
    std::vector<std::mutex> bucket_mutexes(NUM_BUCKETS);
    std::vector<std::vector<size_t>> buckets(NUM_BUCKETS);
    
    for (auto& bucket : buckets) {
        bucket.reserve(total / NUM_BUCKETS); // Pre-allocate approximately
    }
    
    double range_inv = 1.0 / (max_val - min_val);
    
    // Parallel bucket assignment
    const size_t num_threads = get_optimal_thread_count(total);
    std::vector<std::thread> threads;
    
    const size_t chunk_size = total / num_threads;
    
    for (size_t t = 0; t < num_threads; ++t) {
        size_t start = t * chunk_size;
        size_t end = (t == num_threads - 1) ? total : (t + 1) * chunk_size;
        
        threads.emplace_back([start, end, &data, min_val, range_inv, &buckets, &bucket_mutexes]() {
            for (size_t i = start; i < end; ++i) {
                double norm = static_cast<double>(data[i] - min_val) * range_inv;
                int bucket = std::min(NUM_BUCKETS-1, static_cast<int>(norm * NUM_BUCKETS));
                
                // Lock just this bucket for thread safety
                std::lock_guard<std::mutex> lock(bucket_mutexes[bucket]);
                buckets[bucket].push_back(i);
            }
        });
    }
    
    // Wait for all threads to complete
    for (auto& t : threads) {
        if (t.joinable()) t.join();
    }
    threads.clear();
    
    // Sort each bucket in parallel
    for (size_t b = 0; b < NUM_BUCKETS; ++b) {
        if (buckets[b].size() > 32) { // Only sort larger buckets
            threads.emplace_back([&data, b, &buckets]() {
                std::sort(buckets[b].begin(), buckets[b].end(), 
                    [&data](size_t a, size_t b) { return data[a] < data[b]; });
            });
        }
    }
    
    // Wait for all bucket sorts to complete
    for (auto& t : threads) {
        if (t.joinable()) t.join();
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
    
    // Pre-compute ideal positions in parallel
    const size_t num_threads = get_optimal_thread_count(total_elements);
    std::vector<std::thread> threads;
    const size_t chunk_size = total_elements / num_threads;
    
    for (size_t t = 0; t < num_threads; ++t) {
        size_t start = t * chunk_size;
        size_t end = (t == num_threads - 1) ? total_elements : (t + 1) * chunk_size;
        
        // FIXED: Added &dims to the lambda capture
        threads.emplace_back([start, end, &data, min_val, range_inv, &ideals, dims_size, &dims]() {
            for (size_t i = start; i < end; ++i) {
                double norm = static_cast<double>(data[i] - min_val) * range_inv;
                for (size_t d = 0; d < dims_size; ++d) {
                    ideals[i][d] = static_cast<size_t>(norm * (dims[d] - 1));
                }
            }
        });
    }
    
    // Wait for all threads to complete
    for (auto& t : threads) {
        if (t.joinable()) t.join();
    }
    threads.clear();

    // Pre-calculate stride for faster coordinate->index conversion
    std::vector<size_t> strides(dims_size);
    strides[dims_size - 1] = 1;
    for (size_t i = dims_size - 1; i > 0; --i) {
        strides[i - 1] = strides[i] * dims[i];
    }
    
    // Reduced iterations needed due to better starting state
    max_iters = std::max(max_iters / 2, size_t(10)); 
    
    for (size_t iter = 0; iter < max_iters; ++iter) {
        std::atomic<bool> moved{false};
        std::vector<bool> swapped(total_elements, false);
        std::mutex queue_mutex;
        std::priority_queue<std::pair<int, std::pair<size_t, size_t>>> swap_queue;
        
        // Process chunks in parallel to find potential swaps
        const size_t offset = iter % 2;
        std::vector<std::future<void>> futures;
        
        for (size_t t = 0; t < num_threads; ++t) {
            size_t chunk_start = offset + (t * chunk_size * 2);
            size_t chunk_end = std::min(total_elements, offset + ((t+1) * chunk_size * 2));
            
            futures.push_back(std::async(std::launch::async, 
                [chunk_start, chunk_end, &data, &dims, dims_size, &ideals, &swapped, &swap_queue, &queue_mutex, &strides]() {
                    std::vector<size_t> coord(dims_size);
                    std::vector<size_t> new_coord(dims_size);
                    std::vector<std::pair<int, std::pair<size_t, size_t>>> local_swaps;
                    local_swaps.reserve(1000); // Pre-allocate space for local swaps
                    
                    for (size_t idx = chunk_start; idx < chunk_end; idx += 2) {
                        if (idx >= data.size() || swapped[idx]) continue; // FIXED: Changed from dims.size() to data.size()
                        
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

                                    if (neighbor_idx >= data.size() || swapped[neighbor_idx]) continue;

                                    int current_energy = manhattan_distance(coord, ideals[idx]) +
                                                       manhattan_distance(new_coord, ideals[neighbor_idx]);
                                    
                                    int new_energy = manhattan_distance(new_coord, ideals[idx]) +
                                                  manhattan_distance(coord, ideals[neighbor_idx]);

                                    int energy_reduction = current_energy - new_energy;
                                    if (energy_reduction > 1) { // Only consider significant improvements
                                        local_swaps.emplace_back(energy_reduction, std::make_pair(idx, neighbor_idx));
                                    }
                                }
                            }
                        }
                    }
                    
                    // Add local swaps to the global queue at once to reduce lock contention
                    if (!local_swaps.empty()) {
                        std::lock_guard<std::mutex> lock(queue_mutex);
                        for (const auto& swap_item : local_swaps) {
                            swap_queue.push(swap_item);
                        }
                    }
                }
            ));
        }
        
        // Wait for all swap-finding threads
        for (auto& future : futures) {
            future.wait();
        }
        
        // Process swaps (this part is hard to parallelize due to interdependencies)
        // Process only top 80% of swaps for speed
        size_t max_swaps = swap_queue.size() * 0.8;
        size_t swaps_done = 0;
        std::mutex swap_mutex;
        
        while (!swap_queue.empty() && swaps_done < max_swaps) {
            auto [energy_reduction, swap_pair] = swap_queue.top();
            swap_queue.pop();
            
            auto [i, j] = swap_pair;
            
            std::lock_guard<std::mutex> lock(swap_mutex);
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
    
    // Parallel flattening for large matrices
    if (rows * cols > 10000) {
        std::vector<std::vector<T>> thread_results(std::thread::hardware_concurrency());
        std::vector<std::thread> threads;
        
        size_t chunk_size = rows / std::thread::hardware_concurrency();
        if (chunk_size == 0) chunk_size = 1;
        
        for (size_t t = 0; t < std::thread::hardware_concurrency(); ++t) {
            size_t start_row = t * chunk_size;
            size_t end_row = (t == std::thread::hardware_concurrency() - 1) ? rows : (t + 1) * chunk_size;
            
            if (start_row >= rows) break;
            
            threads.emplace_back([start_row, end_row, &matrix, &thread_results, t, cols]() {
                thread_results[t].reserve((end_row - start_row) * cols);
                for (size_t r = start_row; r < end_row; ++r) {
                    thread_results[t].insert(thread_results[t].end(), matrix[r].begin(), matrix[r].end());
                }
            });
        }
        
        for (auto& thread : threads) {
            if (thread.joinable()) thread.join();
        }
        
        for (auto& result : thread_results) {
            flat.insert(flat.end(), result.begin(), result.end());
        }
    } else {
        for (const auto& row : matrix) {
            flat.insert(flat.end(), row.begin(), row.end());
        }
    }
    
    field_space_sort_nd(flat, {rows, cols}, max_iters);
    
    // Parallel reshaping for large matrices
    if (rows * cols > 10000) {
        std::vector<std::thread> threads;
        size_t chunk_size = rows / std::thread::hardware_concurrency();
        if (chunk_size == 0) chunk_size = 1;
        
        for (size_t t = 0; t < std::thread::hardware_concurrency(); ++t) {
            size_t start_row = t * chunk_size;
            size_t end_row = (t == std::thread::hardware_concurrency() - 1) ? rows : (t + 1) * chunk_size;
            
            if (start_row >= rows) break;
            
            threads.emplace_back([start_row, end_row, &matrix, &flat, cols]() {
                for (size_t r = start_row; r < end_row; ++r) {
                    for (size_t c = 0; c < cols; ++c) {
                        matrix[r][c] = flat[r * cols + c];
                    }
                }
            });
        }
        
        for (auto& thread : threads) {
            if (thread.joinable()) thread.join();
        }
    } else {
        size_t idx = 0;
        for (auto& row : matrix) {
            for (auto& val : row) {
                val = flat[idx++];
            }
        }
    }
}

// 3D specialization with parallelism
template<typename T>
void field_space_sort(std::vector<std::vector<std::vector<T>>>& cube, size_t max_iters = 100) {
    if (cube.empty() || cube[0].empty()) return;
    
    const size_t depth = cube.size();
    const size_t rows = cube[0].size();
    const size_t cols = cube[0][0].size();
    const size_t total_elements = depth * rows * cols;
    std::vector<T> flat;
    flat.reserve(total_elements);
    
    // Parallel flattening for large 3D structures
    if (total_elements > 10000) {
        const size_t hw_threads = std::thread::hardware_concurrency();
        std::vector<std::vector<T>> thread_results(hw_threads);
        std::vector<std::thread> threads;
        
        size_t chunk_size = depth / hw_threads;
        if (chunk_size == 0) chunk_size = 1;
        
        for (size_t t = 0; t < hw_threads; ++t) {
            size_t start_depth = t * chunk_size;
            size_t end_depth = (t == hw_threads - 1) ? depth : (t + 1) * chunk_size;
            
            if (start_depth >= depth) break;
            
            threads.emplace_back([start_depth, end_depth, &cube, &thread_results, t, rows, cols]() {
                thread_results[t].reserve((end_depth - start_depth) * rows * cols);
                for (size_t d = start_depth; d < end_depth; ++d) {
                    for (const auto& row : cube[d]) {
                        thread_results[t].insert(thread_results[t].end(), row.begin(), row.end());
                    }
                }
            });
        }
        
        for (auto& thread : threads) {
            if (thread.joinable()) thread.join();
        }
        
        for (auto& result : thread_results) {
            flat.insert(flat.end(), result.begin(), result.end());
        }
    } else {
        for (const auto& layer : cube) {
            for (const auto& row : layer) {
                flat.insert(flat.end(), row.begin(), row.end());
            }
        }
    }
    
    field_space_sort_nd(flat, {depth, rows, cols}, max_iters);
    
    // Parallel reshaping for large 3D structures
    if (total_elements > 10000) {
        const size_t hw_threads = std::thread::hardware_concurrency();
        std::vector<std::thread> threads;
        size_t chunk_size = depth / hw_threads;
        if (chunk_size == 0) chunk_size = 1;
        
        for (size_t t = 0; t < hw_threads; ++t) {
            size_t start_depth = t * chunk_size;
            size_t end_depth = (t == hw_threads - 1) ? depth : (t + 1) * chunk_size;
            
            if (start_depth >= depth) break;
            
            threads.emplace_back([start_depth, end_depth, &cube, &flat, rows, cols]() {
                for (size_t d = start_depth; d < end_depth; ++d) {
                    for (size_t r = 0; r < rows; ++r) {
                        for (size_t c = 0; c < cols; ++c) {
                            cube[d][r][c] = flat[d * rows * cols + r * cols + c];
                        }
                    }
                }
            });
        }
        
        for (auto& thread : threads) {
            if (thread.joinable()) thread.join();
        }
    } else {
        size_t idx = 0;
        for (auto& layer : cube) {
            for (auto& row : layer) {
                for (auto& val : row) {
                    val = flat[idx++];
                }
            }
        }
    }
}

// 4D specialization with parallelism
template<typename T>
void field_space_sort(std::vector<std::vector<std::vector<std::vector<T>>>>& tesseract, size_t max_iters = 100) {
    if (tesseract.empty() || tesseract[0].empty() || tesseract[0][0].empty()) return;
    
    const size_t w = tesseract.size();              // First dimension
    const size_t x = tesseract[0].size();           // Second dimension
    const size_t y = tesseract[0][0].size();        // Third dimension
    const size_t z = tesseract[0][0][0].size();     // Fourth dimension
    const size_t total_elements = w * x * y * z;
    
    std::vector<T> flat;
    flat.reserve(total_elements);
    
    // Parallel flattening for large 4D structures
    if (total_elements > 10000) {
        const size_t hw_threads = std::thread::hardware_concurrency();
        std::vector<std::vector<T>> thread_results(hw_threads);
        std::vector<std::thread> threads;
        
        size_t chunk_size = w / hw_threads;
        if (chunk_size == 0) chunk_size = 1;
        
        for (size_t t = 0; t < hw_threads; ++t) {
            size_t start_w = t * chunk_size;
            size_t end_w = (t == hw_threads - 1) ? w : (t + 1) * chunk_size;
            
            if (start_w >= w) break;
            
            threads.emplace_back([start_w, end_w, &tesseract, &thread_results, t, x, y, z]() {
                thread_results[t].reserve((end_w - start_w) * x * y * z);
                for (size_t i = start_w; i < end_w; ++i) {
                    for (const auto& space : tesseract[i]) {
                        for (const auto& plane : space) {
                            thread_results[t].insert(thread_results[t].end(), plane.begin(), plane.end());
                        }
                    }
                }
            });
        }
        
        for (auto& thread : threads) {
            if (thread.joinable()) thread.join();
        }
        
        for (auto& result : thread_results) {
            flat.insert(flat.end(), result.begin(), result.end());
        }
    } else {
        for (const auto& space : tesseract) {
            for (const auto& plane : space) {
                for (const auto& row : plane) {
                    flat.insert(flat.end(), row.begin(), row.end());
                }
            }
        }
    }
    
    // Sort the flattened data
    field_space_sort_nd(flat, {w, x, y, z}, max_iters);
    
    // Parallel reshaping for large 4D structures
    if (total_elements > 10000) {
        const size_t hw_threads = std::thread::hardware_concurrency();
        std::vector<std::thread> threads;
        size_t chunk_size = w / hw_threads;
        if (chunk_size == 0) chunk_size = 1;
        
        for (size_t t = 0; t < hw_threads; ++t) {
            size_t start_w = t * chunk_size;
            size_t end_w = (t == hw_threads - 1) ? w : (t + 1) * chunk_size;
            
            if (start_w >= w) break;
            
            threads.emplace_back([start_w, end_w, &tesseract, &flat, x, y, z]() {
                for (size_t i = start_w; i < end_w; ++i) {
                    for (size_t j = 0; j < x; ++j) {
                        for (size_t k = 0; k < y; ++k) {
                            for (size_t l = 0; l < z; ++l) {
                                tesseract[i][j][k][l] = flat[i * x * y * z + j * y * z + k * z + l];
                            }
                        }
                    }
                }
            });
        }
        
        for (auto& thread : threads) {
            if (thread.joinable()) thread.join();
        }
    } else {
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
}

#endif // FIELD_SPACE_SORT_HPP