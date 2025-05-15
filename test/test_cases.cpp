#include <cassert>
#include <vector>
#include <iostream>
#include "../include/field_space_sort.hpp"

void test_1d_sort() {
    std::vector<int> arr = {9, 5, 3, 7, 1};
    field_space_sort(arr);
    assert(std::is_sorted(arr.begin(), arr.end()));
}

void test_2d_sort() {
    std::vector<std::vector<int>> matrix = {
        {9, 5}, 
        {3, 1}
    };
    field_space_sort(matrix);
    
    std::vector<int> flat;
    for (const auto& row : matrix) {
        flat.insert(flat.end(), row.begin(), row.end());
    }
    assert(std::is_sorted(flat.begin(), flat.end()));
}

void test_stability() {
    std::vector<float> arr = {1.1, 1.2, 2.1, 2.2};
    field_space_sort(arr);
    assert(arr[0] <= arr[1] && arr[1] <= arr[2] && arr[2] <= arr[3]);
}

int main() {
    test_1d_sort();
    test_2d_sort();
    test_stability();
    std::cout<< "All tests passed!\n";
    return 0;
}