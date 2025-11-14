#ifndef UNIT_OPERATION_BASE_HPP
#define UNIT_OPERATION_BASE_HPP

#include <sundials/sundials_types.h>

#include <cassert>

#include "Components/ComponentSystem.hpp"
#include "EigenDataTypes.hpp"
#include "Reactions/ReactionSystem.hpp"

class Solver;

/// @brief Phase enumeration for multi-phase systems
enum class Phase { Liquid = 0, SolidOrSlurry = 1 };

/// @brief Unit operation types in the process
enum class UnitOperationType { Unknown, Inlet, Outlet, Volume, Pipe, RS_MagneticCaptureProcessChamber };

/**
 * @brief Abstract base for all unit operations
 *
 * Defines common interface (size, indexing, volume access, rhs assembly hook)
 * for concrete operations (Inlet, Outlet, Pipe, Volume, magnetic chamber).
 * Provides generic index mapping for (cell, component, phase) -> y-position.
 */
class UnitOperationBase {
   public:
    UnitOperationBase(const ReactionSystem& reactionSystem, sunindextype n_cells)
        : reactionSystem(reactionSystem), n_cells(n_cells){};

    virtual ~UnitOperationBase() = default;

    virtual UnitOperationType getType() const = 0;
    sunindextype n_components() const { return reactionSystem.componentSystem.n_components; }
    sunindextype y_size() const { return n_cells * n_components(); }

    // for many unit operations, the volume of the liquid phase is constant over time and y, one can call them with default arguments
    virtual const ColVector& get_V_l(realtype t = -1, const realtype* y = nullptr) const = 0;
    virtual realtype errorFunction(realtype t, const realtype* y) const {
        return 0.0;  // Default implementation for unit operations without reactions (Inlet, Outlet)
    };

    virtual sunindextype idx(sunindextype finiteVolumeIdx = 0,
                             sunindextype componentIdx = 0,
                             Phase phase = Phase::Liquid) const {
        return idx<Phase::Liquid>(finiteVolumeIdx, componentIdx, phase);
    };

    // inletIdx and outletIdx are the first indices of the finite volumes of the predecessor and successor unit operations
    // so e.g. [inletIdx, ..., inletIdx + n_components - 1] are the indices of the components in the last finite volume of the predecessor unit operation
    virtual void rhs(realtype t,
                     const realtype* y,
                     realtype* dy_dt,
                     const Solver& solver,
                     bool enable_reactions = true,
                     bool enable_convection = true,
                     bool enable_dispersion = true,
                     bool enable_otherPhysics = true) {
        // Default implementation does nothing, should be overridden in derived classes
    };

    const ReactionSystem& reactionSystem;
    const sunindextype n_cells;  // Number of cells (= finite volumes) in unit operation
    Array y;  // The y vector, used to store the state of the unit operation e.g. for initial conditions, caching, etc. It may be default initialized to size() = 0

   protected:
    // Default function to calculate the index of a component in the y vector
    // used in implementations of addTimeSection methods in derived classes
    // Not virtual, because subclasses with different implementation should use different function signatures
    template <Phase... phases>
    sunindextype idx(sunindextype finiteVolumeIdx = 0, sunindextype componentIdx = 0, Phase phase = Phase::Liquid) const {
        const auto& firstPhaseIdx = finiteVolumeIdx * (sizeof...(phases) * n_components());

        constexpr Phase phaseArray[] = {phases...};
        std::size_t phaseIdx = 0;
        for (; phaseIdx < sizeof...(phases); ++phaseIdx) {
            if (phaseArray[phaseIdx] == phase) break;
        }
        assert(phaseIdx < sizeof...(phases));

        return firstPhaseIdx + phaseIdx * n_components() + componentIdx;
    }

    // y = [y_uOp1, y_uOp2, y_uOp3]
    // y_uOp1 = [y_fV1, y_fV2, y_fV3, y_extra_per_uOp1, y_extra_per_uOp2]
    // y_fV1 = [y_p1, y_p2]
    // y_p1 = [y_c1, y_c2, y_c3]

    // if all equal:
    // n_y = n_unitOp * (n_fv * (n_c*n_p + n_extra_per_fv) + n_extra_per_uOp1)
};

#endif  // UNIT_OPERATION_BASE_HPP