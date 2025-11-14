#ifndef PIPE_HPP
#define PIPE_HPP

#include <functional>

#include "Eigen/Core"
#include "EigenDataTypes.hpp"
#include "UnitOperationBase.hpp"

/**
 * @brief Uniform single-phase pipe unit operation
 *
 * Discretizes a 1D pipe into equally sized finite volume cells (liquid phase only).
 * Provides local RHS contributions for convection (upwind/WENO), axial dispersion,
 * and chemical reactions. Offers ArrayMapper helpers for inlet, outlet and all cells.
 */
class Pipe : public UnitOperationBase {
   public:
    Pipe(const ReactionSystem& reactionSystem,
         sunindextype n_cells,
         realtype crossSectionArea,
         realtype length,
         std::function<realtype(realtype)> flowRateFunction,
         realtype dispersion_coefficient)
        : UnitOperationBase(reactionSystem, n_cells),
          crossSectionArea(crossSectionArea),
          length(length),
          flowRateFunction(flowRateFunction),
          dispersion_coefficient(dispersion_coefficient),
          cell_volume((crossSectionArea * length) / n_cells),
          V_l(ColVector::Constant(n_cells, cell_volume)),
          c(Array(n_cells, reactionSystem.componentSystem.n_components)),
          activities(Array(n_cells, reactionSystem.componentSystem.n_components)),
          dc_dt(Array(n_cells, reactionSystem.componentSystem.n_components)),
          cell_distance(length / n_cells) {}

    UnitOperationType getType() const override { return UnitOperationType::Pipe; }

    void setConstInitialConcentration(const RowVector& concentrations) {
        assert(concentrations.size() == n_components() &&
               "Initial concentration vector size does not match number of components in the component system.");
        y = Array::Zero(n_cells, n_components());
        for (sunindextype i = 0; i < n_cells; ++i) {
            y.row(i) = cell_volume * concentrations;
        }
    };

    const ColVector& get_V_l(realtype t, const realtype* y) const override { return V_l; };

    const ArrayMapper in() const { return ArrayMapper(this, 1, n_components(), 0); }

    const ArrayMapper out() const { return ArrayMapper(this, 1, n_components(), n_components() * (n_cells - 1)); }

    const ArrayMapper all() const { return ArrayMapper(this, n_cells, n_components(), 0); }

    realtype errorFunction(realtype t, const realtype* y) const override;

    void rhs(realtype t,
             const realtype* y,
             realtype* dy_dt,
             const Solver& solver,
             bool enable_reactions = true,
             bool enable_convection = true,
             bool enable_dispersion = true,
             bool enable_otherPhysics = true) override;

    const realtype crossSectionArea;                           // Cross-sectional area of the pipe
    const realtype length;                                     // Length of the pipe
    const std::function<realtype(realtype)> flowRateFunction;  // returns flow rate for a given time
    const realtype dispersion_coefficient;                     // Dispersion coefficient for the pipe

    // Derived helper parameters
    const realtype cell_distance;  // Distance between the centers of two adjacent cells
    const realtype cell_volume;    // Volume of each finite volume in the pipe

   private:
    const ColVector V_l;
    mutable Array c;
    mutable Array activities;
    mutable Array dc_dt;
};

#endif  // PIPE_HPP