//
// Created by cui on 25-6-6.
//
#ifndef TOOLS_H
#define TOOLS_H
#include <chrono>
#include<iostream>
#include <algorithm>
#include <cstring>
#include <fstream>

#include <vector>
#include <ctime>
#include <string>

#define FF 4
using namespace std;




int ctz128_uint128(__uint128_t num) {
    uint64_t low_bits = num;
    uint64_t high_bits = num >> 64;

    if (low_bits != 0) {
        return __builtin_ctzll(low_bits);
    } else if (high_bits != 0) {
        return 64 + __builtin_ctzll(high_bits);
    } else {
        return 128; // 如果输入为 0，则返回 128
    }
}

int clz128_uint128(__uint128_t num) {
    uint64_t low_bits = num;
    uint64_t high_bits = num >> 64;

    if (high_bits != 0) {
        return __builtin_clzll(high_bits);
    } else if (low_bits != 0) {
        return 64 + __builtin_clzll(low_bits);
    } else {
        return 128; // 如果输入为 0，则返回 128
    }
}

// 为 unsigned __int128 类型定义 operator<<
std::ostream &operator<<(std::ostream &os, const unsigned __int128 &value) {
    // 由于 __int128 是一个非常大的整数，我们需要转换成字符串来输出
    std::string output;
    unsigned __int128 temp = value;
    do {
        output += '0' + (temp % 10);
        temp /= 10;
    } while (temp != 0);
    std::reverse(output.begin(), output.end());
    os << output;
    return os;
}

int countSetBits128(const __uint128_t &n) {
    uint64_t low_part = static_cast<uint64_t>(n); // 低位部分
    uint64_t high_part = static_cast<uint64_t>(n >> 64); // 高位部分
    return __builtin_popcountll(low_part) + __builtin_popcountll(high_part);
}
#endif //TOOLS_H
