#include "UnitOperations/Volume.hpp"

#include "ConvectionDispersionOperators.hpp"
#include "Eigen/Core"
#include "Solver.hpp"
#include "UnitOperations/Inlet.hpp"

realtype Volume::errorFunction(realtype t, const realtype* y) const {
    if (std::isnan(cell_volume)) {
        return 0.0;  // No reactions for open reservoirs
    }

    PhaseStride stride(n_components(), 1);
    ConstArrayMap m(y, n_cells, n_components(), stride);

    c = (1 / cell_volume) * m;
    c = c.cwiseMax(0);  // Clamp negatives to 0

    activities.setZero();

    return reactionSystem.errorFunction(t, c, activities);
}

void Volume::rhs(realtype t,
                 const realtype* y,
                 realtype* dy_dt,
                 const Solver& solver,
                 bool enable_reactions,
                 bool enable_convection,
                 bool enable_dispersion,
                 bool enable_otherPhysics) {
    if (std::isnan(cell_volume)) {
        return;  // No internal dynamics for open reservoirs
    }

    PhaseStride stride(n_components(), 1);

    ConstArrayMap m(y, n_cells, n_components(), stride);
    ArrayMap dm_dt(dy_dt, n_cells, n_components(), stride);

    if (enable_reactions) {
        c = (1 / cell_volume) * m;
#if LOG_ENABLED
        if ((V_l.array() <= 0).any()) {
            std::ostringstream oss;
            oss << "Warning: Negative value detected in V_l in Pipe at time t = " << t << ". V_l = " << V_l.transpose()
                << "\n";
            LOG("unallowed_values.log", oss.str());
        }
#endif
        c = c.cwiseMax(0);  // Clamp negatives to 0
        // c = c.cwiseAbs();     // Take absolute values

        activities.setZero();
        dc_dt.setZero();

        reactionSystem.rhs(t, c, activities, dc_dt);

        dm_dt += cell_volume * dc_dt;
    }
}
