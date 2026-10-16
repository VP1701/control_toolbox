// matrix_io.hpp
// Reads matrices written by matlab_testcode/export_ss.m (host-side tests only).
//
// File format (whitespace separated, '#' starts a comment):
//     <name> <rows> <cols>
//     <rows * cols values, row-major>
//     ...
// Values are written with %.17g, so doubles round-trip exactly.

#ifndef IMPORT_SS_HPP
#define IMPORT_SS_HPP

#include <cerrno>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "matrix.h"
#include "status.hpp"

namespace matrix_io {

using MatrixMap = std::map<std::string, matrix>;

namespace detail {

// Upper limit for one dimension, rejects corrupt headers before allocating.
constexpr long MAX_DIM = 100000;

inline std::vector<std::string> tokenize(std::ifstream& file) {
    std::vector<std::string> tokens;
    std::string line;
    while (std::getline(file, line)) {
        const std::size_t comment = line.find('#');
        if (comment != std::string::npos) {
            line.erase(comment);
        }
        std::istringstream stream(line);
        std::string token;
        while (stream >> token) {
            tokens.push_back(token);
        }
    }
    return tokens;
}

inline bool is_name(const std::string& token) {
    // a name starts with a letter or '_', which a number never does
    // (except inf/nan, which are rejected as names by this rule too)
    const char first = token[0];
    const bool letter = (first >= 'A' && first <= 'Z') || (first >= 'a' && first <= 'z');
    const bool numeric_word = (token == "Inf" || token == "inf" ||
                               token == "NaN" || token == "nan");
    return (letter || first == '_') && !numeric_word;
}

inline bool to_dim(const std::string& token, int& value) {
    char* end = nullptr;
    errno = 0;
    const long parsed = std::strtol(token.c_str(), &end, 10);
    if (end == token.c_str() || *end != '\0' || errno == ERANGE ||
        parsed < 0 || parsed > MAX_DIM) {
        return false;
    }
    value = static_cast<int>(parsed);
    return true;
}

inline bool to_double(const std::string& token, double& value) {
    // strtod also accepts Inf, -Inf and NaN as written by MATLAB
    char* end = nullptr;
    errno = 0;
    value = std::strtod(token.c_str(), &end);
    const bool consumed_all = (end != token.c_str()) && (*end == '\0');
    const bool overflow = (errno == ERANGE) && std::isinf(value);
    return consumed_all && !overflow;
}

}  // namespace detail

// Reads every matrix in the file. 'out' is only modified on success.
inline Status load_matrices(const std::string& filename, MatrixMap& out) {
    std::ifstream file(filename);
    if (!file) {
        return Status::FileError;
    }

    const std::vector<std::string> tokens = detail::tokenize(file);
    MatrixMap result;
    std::size_t pos = 0;

    while (pos < tokens.size()) {
        if (tokens.size() - pos < 3) {
            return Status::ParseError;
        }
        const std::string& name = tokens[pos];
        int rows = 0;
        int cols = 0;
        if (!detail::is_name(name) || !detail::to_dim(tokens[pos + 1], rows) ||
            !detail::to_dim(tokens[pos + 2], cols) || result.count(name) != 0) {
            return Status::ParseError;
        }
        pos += 3;

        const std::size_t count = static_cast<std::size_t>(rows) * static_cast<std::size_t>(cols);
        if (count > tokens.size() - pos) {
            return Status::ParseError;
        }

        matrix m(rows, cols);
        for (std::size_t i = 0; i < count; ++i) {
            if (!detail::to_double(tokens[pos], m.data[i])) {
                return Status::ParseError;
            }
            ++pos;
        }
        result[name] = m;
    }

    out = result;
    return Status::Ok;
}

// Copies one matrix out of the map.
inline Status get(const MatrixMap& matrices, const std::string& name, matrix& out) {
    const auto it = matrices.find(name);
    if (it == matrices.end()) {
        return Status::NotFound;
    }
    out = it->second;
    return Status::Ok;
}

// Reads A, B, C, D and checks that the sizes form a valid state-space model.
// The outputs are only modified on success.
inline Status load_ss(const std::string& filename,
                      matrix& A, matrix& B, matrix& C, matrix& D) {
    MatrixMap matrices;
    Status s = load_matrices(filename, matrices);
    if (s != Status::Ok) {
        return s;
    }

    matrix a;
    matrix b;
    matrix c;
    matrix d;
    if (get(matrices, "A", a) != Status::Ok || get(matrices, "B", b) != Status::Ok ||
        get(matrices, "C", c) != Status::Ok || get(matrices, "D", d) != Status::Ok) {
        return Status::NotFound;
    }

    const int n = a.rows;
    const int m = b.columns;
    const int p = c.rows;
    if (a.columns != n || b.rows != n || c.columns != n || d.rows != p || d.columns != m) {
        return Status::DimensionMismatch;
    }

    A = a;
    B = b;
    C = c;
    D = d;
    return Status::Ok;
}

}  // namespace import_ss

#endif  // IMPORT_SS_HPP