#ifndef VOLUME_HPP
#define VOLUME_HPP

#include <sundials/sundials_types.h>

#include <functional>

#include "Eigen/Core"
#include "EigenDataTypes.hpp"
#include "UnitOperationBase.hpp"

/**
 * @brief Single-cell volume unit operation
 *
 * Models a homogeneous control volume (one finite volume cell) with optional
 * fixed liquid volume. Supports reactions; useful for dead volumes, reservoirs,
 * or fraction collection.
 */
class Volume : public UnitOperationBase {
   public:
    Volume(const ReactionSystem& reactionSystem, realtype volume = std::numeric_limits<realtype>::quiet_NaN())
        : UnitOperationBase(reactionSystem, 1),
          cell_volume(volume),
          V_l(ColVector::Constant(n_cells, cell_volume)),
          c(Array(n_cells, n_components())),
          activities(Array(n_cells, n_components())),
          dc_dt(Array(n_cells, n_components())) {}

    UnitOperationType getType() const override { return UnitOperationType::Volume; }

    void setConstInitialConcentration(const RowVector& concentrations) {
        if (V_l.hasNaN()) {
            throw std::runtime_error("Cell volume is not set for Unit Operation 'Volume'.");
        } else {
            assert(concentrations.size() == n_components() &&
                   "Initial concentration vector size does not match number of components in the component system.");
            y = Array::Zero(n_cells, n_components());
            for (sunindextype i = 0; i < n_cells; ++i) {
                y.row(i) = cell_volume * concentrations;
            }
        }
    };

    // Volume suppports both fixed volumes e.g. dead volumes and volumes without a precise size e.g. open reservoirs
    virtual const ColVector& get_V_l(realtype t, const realtype* y) const override {
        if (V_l.hasNaN()) {
            throw std::runtime_error("Cell volume is not set for Unit Operation 'Volume'.");
        } else {
            return V_l;
        }
    };

    const ArrayMapper in() const { return ArrayMapper(this, 1, n_components(), 0); }
    const ArrayMapper out() const { return ArrayMapper(this, 1, n_components(), 0); }
    const ArrayMapper all() const { return ArrayMapper(this, 1, n_components(), 0); }

    realtype errorFunction(realtype t, const realtype* y) const override;

    void rhs(realtype t,
             const realtype* y,
             realtype* dy_dt,
             const Solver& solver,
             bool enable_reactions = true,
             bool enable_convection = true,
             bool enable_dispersion = true,
             bool enable_otherPhysics = true) override;

    const realtype cell_volume;  // Volume of each finite volume in the pipe

   private:
    const ColVector V_l;
    mutable Array c;
    mutable Array activities;
    mutable Array dc_dt;
};

#endif  // VOLUME_HPP