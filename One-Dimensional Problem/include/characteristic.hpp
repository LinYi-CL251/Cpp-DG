#pragma once

#include <array>
#include <tuple>

namespace dg {

template<class State>
struct CharacteristicBasis {
    using Matrix = std::array<State, std::tuple_size<State>::value>;
    Matrix left{};  // Rows are left eigenvectors.
    Matrix right{}; // Columns are right eigenvectors.

    static CharacteristicBasis identity() {
        CharacteristicBasis basis;
        for (std::size_t j = 0; j < std::tuple_size<State>::value; ++j) {
            basis.left[j][j] = basis.right[j][j] = 1.0;
        }
        return basis;
    }
};

template<class State>
State matrix_vector(const typename CharacteristicBasis<State>::Matrix& matrix,
                    const State& vector) {
    State result{};
    for (std::size_t row = 0; row < vector.size(); ++row) {
        for (std::size_t col = 0; col < vector.size(); ++col) {
            result[row] += matrix[row][col] * vector[col];
        }
    }
    return result;
}

} // namespace dg