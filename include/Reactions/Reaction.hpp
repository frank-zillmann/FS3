#ifndef REACTION_HPP
#define REACTION_HPP

#include <functional>
#include <memory>

#include "EigenDataTypes.hpp"
#include "sundials/sundials_types.h"

/**
 * @brief Vectorized reaction wrapper
 *
 * Holds RHS functor operating on (concentrations, activities) over multiple cells
 * and optional error diagnostic (e.g. equilibrium deviation or pH). Constructed
 * via factory helpers for mass action variants.
 */
struct Reaction {
    // Eigen-based interface for vectorized operations across multiple cells
    Reaction(
        std::function<void(realtype t, const ConstArrayMap& concentrations, const ConstArrayMap& activities, ArrayMap& dc_dt)> rhsFunction,
        std::function<realtype(realtype t, const ConstArrayMap& concentrations, const ConstArrayMap& activities)> errorFunction =
            [](realtype, const ConstArrayMap&, const ConstArrayMap&) { return 0.0; })
        : rhs(rhsFunction), errorFunction(errorFunction) {}

    std::function<void(realtype t, const ConstArrayMap& concentrations, const ConstArrayMap& activities, ArrayMap& dc_dt)> rhs;

    std::function<realtype(realtype t, const ConstArrayMap& concentrations, const ConstArrayMap& activities)> errorFunction;
};

using ReactionHook = std::function<void(realtype t, ArrayMap& concentrations, ArrayMap& activities, ArrayMap& dc_dt)>;

inline Reaction wrap_reaction(const Reaction& base, ReactionHook pre = {}, ReactionHook after = {}) {
    if (!pre && !after) {
        return base;
    }

    auto concentrations_buffer = std::make_shared<Array>();
    auto activities_buffer = std::make_shared<Array>();
    auto dc_dt_buffer = std::make_shared<Array>();

    return Reaction(
        [base, pre, after, concentrations_buffer, activities_buffer,
         dc_dt_buffer](realtype t, const ConstArrayMap& concentrations, const ConstArrayMap& activities, ArrayMap& dc_dt) {
            const auto rows = activities.rows();
            const auto cols = activities.cols();

            if (concentrations_buffer->rows() != rows || concentrations_buffer->cols() != cols) {
                concentrations_buffer->resize(rows, cols);
            }
            if (activities_buffer->rows() != rows || activities_buffer->cols() != cols) {
                activities_buffer->resize(rows, cols);
            }
            if (dc_dt_buffer->rows() != rows || dc_dt_buffer->cols() != cols) {
                dc_dt_buffer->resize(rows, cols);
            }

            *concentrations_buffer = concentrations;
            *activities_buffer = activities;
            dc_dt_buffer->setZero();

            ArrayMap concentrations_map(concentrations_buffer->data(), rows, cols, PhaseStride(cols, 1));
            ArrayMap activities_map(activities_buffer->data(), rows, cols, PhaseStride(cols, 1));
            ArrayMap dc_dt_map(dc_dt_buffer->data(), rows, cols, PhaseStride(cols, 1));
            ConstArrayMap concentrations_const_map(concentrations_buffer->data(), rows, cols, PhaseStride(cols, 1));
            ConstArrayMap activities_const_map(activities_buffer->data(), rows, cols, PhaseStride(cols, 1));

            if (pre) {
                pre(t, concentrations_map, activities_map, dc_dt_map);
            }

            base.rhs(t, concentrations_const_map, activities_const_map, dc_dt_map);

            if (after) {
                after(t, concentrations_map, activities_map, dc_dt_map);
            }
            dc_dt += dc_dt_map;
        },
        [base, pre, concentrations_buffer, activities_buffer,
         dc_dt_buffer](realtype t, const ConstArrayMap& concentrations, const ConstArrayMap& activities) -> realtype {
            if (!pre) {
                return base.errorFunction(t, concentrations, activities);
            }

            const auto rows = activities.rows();
            const auto cols = activities.cols();

            if (concentrations_buffer->rows() != rows || concentrations_buffer->cols() != cols) {
                concentrations_buffer->resize(rows, cols);
            }
            if (activities_buffer->rows() != rows || activities_buffer->cols() != cols) {
                activities_buffer->resize(rows, cols);
            }
            if (dc_dt_buffer->rows() != rows || dc_dt_buffer->cols() != cols) {
                dc_dt_buffer->resize(rows, cols);
            }

            *concentrations_buffer = concentrations;
            *activities_buffer = activities;
            dc_dt_buffer->setZero();

            ArrayMap concentrations_map(concentrations_buffer->data(), rows, cols, PhaseStride(cols, 1));
            ArrayMap activities_map(activities_buffer->data(), rows, cols, PhaseStride(cols, 1));
            ArrayMap dc_dt_map(dc_dt_buffer->data(), rows, cols, PhaseStride(cols, 1));
            ConstArrayMap concentrations_const_map(concentrations_buffer->data(), rows, cols, PhaseStride(cols, 1));
            ConstArrayMap activities_const_map(activities_buffer->data(), rows, cols, PhaseStride(cols, 1));

            pre(t, concentrations_map, activities_map, dc_dt_map);
            return base.errorFunction(t, concentrations_const_map, activities_const_map);
        });
}

#endif  // REACTION_HPP