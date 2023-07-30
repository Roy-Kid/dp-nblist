#pragma once // 使用#pragma once确保头文件只被编译一次

#include <iostream>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/sequence.h>
#include <thrust/scatter.h>
#include <thrust/unique.h>

// 需要类模板化
template <typename T>
struct Dynamic2DArray {
    thrust::device_vector<int> helper; // 辅助数组，记录每行的值数量
    thrust::device_vector<T> flat_data; // 扁平化的一维数组用于存储二维数组的元素
    int max_rows; // 最大行数
    int max_cols; // 最大列数

    // Default constructor
    Dynamic2DArray() : max_rows(0), max_cols(0) {
        helper.resize(0);
        flat_data.resize(0);
    }

    // 构造函数，初始化数组
    Dynamic2DArray(int rows, int cols) : max_rows(rows), max_cols(cols) {
        helper.resize(rows, 0);
        for(int i = 0; i < rows; i++) {
            helper[i] = 0;
        }

        flat_data.resize(rows * cols, 0);
    }

    // push_back操作，在指定行添加一个元素
    void push_back(int row, T value) {
        if (row < max_rows) {
            flat_data[row*max_cols + helper[row]] = value; // 将新元素添加到扁平化的一维数组中
            helper[row]++; // 更新辅助数组，增加新行的大小
        }
    }

    // 获取二维数组元素
    T getElement(int row, int col) {
        if (row < max_rows && col < max_cols) {
            return flat_data[row*max_cols + col];
        }
        return T{}; // 返回默认值（可以根据需求修改）
    }

    thrust::device_vector<T> getVector(int row) {
        thrust::device_vector<T> res(helper[row]);
        for (int i = 0; i < helper[row]; i++) {
            res[i] = flat_data[i + row * max_cols];
        }
        return res;
    }

    void resize(int row, int col) {
        max_rows = row;      
        max_cols = col;

        helper.resize(row);
        flat_data.resize(row * col);
        // std::cout << "sign" << std::endl;
        // std::cout << max_cols << " " << max_cols << std::endl;
        // std::cout << row << " " << col << std::endl;
    }

    // 输出二维数组内容
    void print() {
        // std::cout << max_rows << " " << max_cols << std::endl;
        for (int i = 0; i < max_rows; i++) {
            int start = i * max_cols;
            int end = helper[i] + start;
            std::cout << "row: " << i << std::endl;
            for (int j = start; j < end; j++) {
                std::cout << flat_data[j] << " ";
            }
            std::cout << std::endl;
        }
    }

    // 输出指定行的二维数组内容
    void print(int rowNum) {
        // std::cout << max_rows << " " << max_cols << std::endl;
        int start = rowNum * max_cols;
        int end = helper[rowNum] + start;
        for (int j = start; j < end; j++) {
            std::cout << flat_data[j] << " ";
        }
        std::cout << std::endl;
    }
};

thrust::device_vector<int> vector2set_int(const thrust::device_vector<int>& d_vec) {
    thrust::device_vector<int> d_vec_copy = d_vec;

    // 去除重复
    // thrust::sort(d_vec_copy.begin(), d_vec_copy.end());
    thrust::device_vector<int>::iterator new_end = thrust::unique(d_vec_copy.begin(), d_vec_copy.end());

    // 重置大小
    d_vec_copy.resize(new_end - d_vec_copy.begin());

    return d_vec_copy;
}




// int main() {
//     int rows = 3;
//     int max_cols = 5;
//     Dynamic2DArray array;
//     array(rows, max_cols);

//     // 进行 push_back 操作，在第 1 行添加一个元素 100
//     array.push_back(1, 31);
//     array.push_back(1, 31);

//     array.push_back(0, 11);
//     array.push_back(0, 12);
//     array.push_back(0, 13);

//     array.push_back(2, 21);
//     array.push_back(2, 22);
//     array.push_back(2, 23);

//     // 输出结果
//     array.print();
//     // array.print(0);
//     // array.print(1);
//     // array.print(2);

//     for(int i = 0; i < rows; i++){
//         for(int j = 0; j < max_cols; j++) {
//             std::cout << array.getElement(i, j) << " ";
//         }
//     }

//     return 0;
// }
