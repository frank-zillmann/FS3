#ifndef REACTION_HPP
#define REACTION_HPP

#include <functional>
#include <optional>

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

#endif  // REACTION_HPP