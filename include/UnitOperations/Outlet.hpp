#ifndef OUTLET_HPP
#define OUTLET_HPP

#include "UnitOperationBase.hpp"

/**
 * @brief Boundary outlet unit operation
 *
 * Collects outflow from connected units. No internal finite volumes;
 * provides an ArrayMapper for its inlet. Can be segmented over time.
 */
class Outlet : public UnitOperationBase {
   public:
    Outlet(const ReactionSystem& reactionSystem) : UnitOperationBase(reactionSystem, 0) {}

    virtual UnitOperationType getType() const override { return UnitOperationType::Outlet; }

    const ColVector& get_V_l(realtype t, const realtype* y) const override { return V_l; }

    void addTimeSection(realtype endTime,
                        const std::vector<UnitOperationBase*>& adjacentUnitOperations,
                        const std::vector<realtype>& flowRates);

    const ArrayMapper in() const { return ArrayMapper(this, 1, n_components(), 0); }

    sunindextype idx(sunindextype finiteVolumeIdx, sunindextype componentIdx, Phase phase) const override {
        throw std::logic_error("Outlet does not have finite volumes, idx should not be called.");
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

    const ColVector V_l = ColVector::Constant(1, 1.0);  // makes c = m
};

#endif  // OUTLET_HPP