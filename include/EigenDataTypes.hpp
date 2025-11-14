#ifndef EIGENDATATYPES_HPP
#define EIGENDATATYPES_HPP

#include <concepts>
#include <cstddef>
#include <functional>

#include "Eigen/Dense"
#include "sundials/sundials_types.h"

// Forward declaration
class UnitOperationBase;

// dynamic-sized Vector
using Vector = Eigen::Vector<realtype, Eigen::Dynamic>;
using VectorMap = Eigen::Map<Vector, 1>;
using ConstVectorMap = Eigen::Map<const Vector, 1>;

// dynamic-sized 2D Array, row‐major
using Array = Eigen::Array<realtype, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

// stride type for (rowStride / outer stride / how many elements between two rows, colStride / inner stride / how many elements between two columns withing a row)
using PhaseStride = Eigen::Stride<Eigen::Dynamic, 1>;
// Strided Maps
using ArrayMap = Eigen::Map<Array, 0, PhaseStride>;
using ConstArrayMap = Eigen::Map<const Array, 0, PhaseStride>;

// Column and row *array* types (elementwise semantics)
using ColVector = Eigen::Array<realtype, Eigen::Dynamic, 1>;
using RowVector = Eigen::Array<realtype, 1, Eigen::Dynamic>;

// Non-strided (contiguous) maps for vectors
using ColVectorMap = Eigen::Map<ColVector>;
using ConstColVectorMap = Eigen::Map<const ColVector>;
using RowVectorMap = Eigen::Map<RowVector>;
using ConstRowVectorMap = Eigen::Map<const RowVector>;

/**
 * @brief Lightweight mapper for strided Array views into y
 *
 * Wraps raw pointers to create Eigen::Array maps with configurable row stride
 * to support phase-interleaved layouts. Used to map unit-operation slices of
 * the global state without copies.
 */
class ArrayMapper {
   public:
    // Default constructor - initializes with sensible defaults
    ArrayMapper() : unitOperation_ptr(nullptr), rows(0), cols(0), startIdx(0), stride(1, 1) {}

    ArrayMapper(const UnitOperationBase* unitOperation_ptr,
                const Eigen::Index rows,
                const Eigen::Index cols,
                const Eigen::Index startIdx = 0,
                const Eigen::Index n_phases = 1)
        : unitOperation_ptr(unitOperation_ptr), rows(rows), cols(cols), startIdx(startIdx), stride(cols * n_phases, 1) {}

    ArrayMapper(const Eigen::Index rows,
                const Eigen::Index cols,
                const Eigen::Index startIdx = 0,
                const Eigen::Index n_phases = 1)
        : unitOperation_ptr(nullptr), rows(rows), cols(cols), startIdx(startIdx), stride(cols * n_phases, 1) {}

    // Pointer to the UnitOperation, nullptr if not associated with one
    const UnitOperationBase* unitOperation_ptr;

    Eigen::Index startIdx;
    Eigen::Index rows, cols;
    PhaseStride stride;

    ArrayMap operator()(realtype* data) const { return ArrayMap(data + startIdx, rows, cols, stride); }
    ConstArrayMap operator()(const realtype* data) const { return ConstArrayMap(data + startIdx, rows, cols, stride); }

    // Avoid usage because capturing the function is slow
    // ArrayMap operator()(realtype* data,
    //                     std::function<std::size_t(const UnitOperationBase*)>& getUnitOperationStartIdx) const {
    //     return ArrayMap(data + getUnitOperationStartIdx(unitOperation_ptr) + startIdx, rows, cols, stride);
    // }
    // ConstArrayMap operator()(const realtype* data,
    //                          std::function<std::size_t(const UnitOperationBase*)>& getUnitOperationStartIdx) const {
    //     return ConstArrayMap(data + getUnitOperationStartIdx(unitOperation_ptr) + startIdx, rows, cols, stride);
    // }
};

// Concepts for better error messages and type constraints
template <typename T>
concept EigenArrayLike = requires(T& t, const T& ct) {
    typename T::Scalar;
    { ct.rows() } -> std::convertible_to<Eigen::Index>;
    { ct.cols() } -> std::convertible_to<Eigen::Index>;
    { t.topRows(1) };
    { t.bottomRows(1) };
    { t.middleRows(1, 1) };
    { ct.topRows(1) };
    { ct.bottomRows(1) };
    { ct.middleRows(1, 1) };
    { t.row(0) };
    { ct.row(0) };
    { t.col(0) };
    { ct.col(0) };
};

template <typename T>
concept WritableEigenArrayLike = EigenArrayLike<T> && requires(T& t) {
    t.row(0) -= t.row(0);
    t.row(0) += t.row(0);
    t.col(0) *= t.col(0);
    t.col(0) /= t.col(0);
};

#endif  // EIGENDATATYPES_HPP