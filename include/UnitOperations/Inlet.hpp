#ifndef INLET_HPP
#define INLET_HPP

#include <sundials/sundials_types.h>

#include <cmath>
#include <functional>
#include <limits>
#include <vector>

#include "EigenDataTypes.hpp"
#include "UnitOperationBase.hpp"

/**
 * @brief Boundary inlet unit operation
 *
 * Supplies a time-dependent boundary composition via a user function.
 * Has no internal finite volumes; provides an ArrayMapper for its outlet.
 */
class Inlet : public UnitOperationBase {
   public:
    Inlet(const ReactionSystem& reactionSystem, std::function<RowVector(const realtype&)> solution_function)
        : UnitOperationBase(reactionSystem, 0), solution_function(solution_function){};

    UnitOperationType getType() const override { return UnitOperationType::Inlet; }

    const ColVector& get_V_l(realtype t, const realtype* y) const override { return V_l; }

    const ArrayMapper out() const { return ArrayMapper(this, 1, n_components(), 0); }

    RowVector getSolution(const realtype& time) const { return solution_function(time); }

    sunindextype idx(sunindextype finiteVolumeIdx, sunindextype componentIdx, Phase phase) const override {
        throw std::logic_error("Inlet does not have finite volumes, idx should not be called.");
        return 0;  // This line will never be reached, but is needed to avoid compiler warnings
    }

    void rhs(realtype t,
             const realtype* y,
             realtype* dy_dt,
             const Solver& solver,
             bool enable_reactions = true,
             bool enable_convection = true,
             bool enable_dispersion = true,
             bool enable_otherPhysics = true) override {
        return;
    }

    const std::function<RowVector(const realtype&)> solution_function;
    const ColVector V_l = ColVector::Constant(1, 1.0);  // makes c = m
};

#endif  // INLET_HPP