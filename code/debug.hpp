#pragma once

#include <array>
#include <iostream>

template <size_t N>
std::ostream& operator<<(std::ostream& os, const std::array<int, N>& a) {
    bool fst = true;
    os << "[";
    for (auto x : a) {
        if (fst)
            fst = false;
        else
            os << ", ";
        os << x;
    }
    os << "]";
    return os;
}

