#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_erkstep.h>
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>
#include <sunlinsol/sunlinsol_band.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_band.h>

#include <algorithm>
#include <iomanip>
#include <ios>
#include <iostream>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#include "Components/Component.hpp"
#include "EigenDataTypes.hpp"
#include "Logger.hpp"
#include "Process.hpp"
#include "Solver.hpp"
#include "cnpy.h"

int Solver::rhs(realtype t, N_Vector y_sundials, N_Vector dy_dt_sundials, void* user_data) {
    // CAREFUL: This function is static, so it cannot access non-static members directly.

    Solver* solver = static_cast<Solver*>(user_data);
    const realtype* y = N_VGetArrayPointer(y_sundials);
    realtype* dy_dt = N_VGetArrayPointer(dy_dt_sundials);

#if LOG_ENABLED
    auto& log_state = solver->rhs_log_state;
    if (!log_state.is_initialized()) {
        log_state.initialize(solver->ySize, solver->process.unitOperations.size());
    }
    log_state.y_copy.assign(y, y + solver->ySize);
#endif

    // Zero-initialize dy_dt before accumulating contributions
    std::fill(dy_dt, dy_dt + solver->ySize, 0.0);

    std::vector<realtype> t_unitOperations;
    t_unitOperations.reserve(solver->process.unitOperations.size());
    for (auto const& unitOperation_ptr : solver->process.unitOperations) {
        size_t first_idx = solver->getUnitOperationStartIdx(unitOperation_ptr.get());
        realtype t_unitOperation = 0.0;
        BENCHMARK(t_unitOperation,
                  { unitOperation_ptr->rhs(t, y + first_idx, dy_dt + first_idx, *solver, true, true, true, true); });
        t_unitOperations.push_back(t_unitOperation);
    }
    realtype t_connections = 0.0;
    BENCHMARK(t_connections, { solver->process.rhs_connections(t, y, dy_dt, *solver); });

#if BENCHMARK_ENABLED
    {
        std::ostringstream oss;
        oss.setf(std::ios::scientific);
        oss.precision(6);
        realtype t_total = std::accumulate(t_unitOperations.begin(), t_unitOperations.end(), 0.0) + t_connections;
        oss << "t=" << t << "\t\tt_total=" << (t_total) << "\t\tt_conncections=" << (t_connections);
        for (size_t i = 0; i < t_unitOperations.size(); ++i) {
            const auto& uo_ptr = solver->process.unitOperations[i];
            oss << "\t\tt_unitOperation[" << i << ",type=" << static_cast<int>(uo_ptr->getType())
                << "]=" << (t_unitOperations[i]);
        }
        oss << "\n";
        LOG_BENCHMARK("rhs.log", oss.str());
    }
#endif

#if LOG_ENABLED
    log_rhs_call(t, y_sundials, dy_dt_sundials, log_state, "rhs.log");
#endif

    return 0;
}

int Solver::rhs_nonStiff(realtype t, N_Vector y_sundials, N_Vector dy_dt_sundials, void* user_data) {
    // CAREFUL: This function is static, so it cannot access non-static members directly.

    Solver* solver = static_cast<Solver*>(user_data);
    const realtype* y = N_VGetArrayPointer(y_sundials);
    realtype* dy_dt = N_VGetArrayPointer(dy_dt_sundials);

#if LOG_ENABLED
    auto& log_state = solver->rhs_nonStiff_log_state;
    if (!log_state.is_initialized()) {
        log_state.initialize(solver->ySize, solver->process.unitOperations.size());
    }
    log_state.y_copy.assign(y, y + solver->ySize);
#endif

    // Zero-initialize dy_dt before accumulating contributions
    std::fill(dy_dt, dy_dt + solver->ySize, 0.0);

    std::vector<realtype> t_unitOperations;
    t_unitOperations.reserve(solver->process.unitOperations.size());
    for (auto const& unitOperation_ptr : solver->process.unitOperations) {
        size_t first_idx = solver->getUnitOperationStartIdx(unitOperation_ptr.get());
        realtype t_unitOperation = 0.0;
        BENCHMARK(t_unitOperation,
                  { unitOperation_ptr->rhs(t, y + first_idx, dy_dt + first_idx, *solver, false, true, true, true); });
        t_unitOperations.push_back(t_unitOperation);
    }
    realtype t_connections = 0.0;
    BENCHMARK(t_connections, { solver->process.rhs_connections(t, y, dy_dt, *solver); });

#if BENCHMARK_ENABLED
    {
        std::ostringstream oss;
        oss.setf(std::ios::scientific);
        oss.precision(6);
        realtype t_total = std::accumulate(t_unitOperations.begin(), t_unitOperations.end(), 0.0) + t_connections;
        oss << "t=" << t << "\t\tt_total=" << (t_total) << "\t\tt_conncections=" << (t_connections);
        for (size_t i = 0; i < t_unitOperations.size(); ++i) {
            const auto& uo_ptr = solver->process.unitOperations[i];
            oss << "\t\tt_unitOperation[" << i << ",type=" << static_cast<int>(uo_ptr->getType())
                << "]=" << (t_unitOperations[i]);
        }
        oss << "\n";
        LOG_BENCHMARK("rhs_nonStiff.log", oss.str());
    }
#endif

#if LOG_ENABLED
    log_rhs_call(t, y_sundials, dy_dt_sundials, log_state, "rhs_nonStiff.log");
#endif

    return 0;
}

int Solver::rhs_stiff(realtype t, N_Vector y_sundials, N_Vector dy_dt_sundials, void* user_data) {
    // CAREFUL: This function is static, so it cannot access non-static members directly.

    Solver* solver = static_cast<Solver*>(user_data);
    const realtype* y = N_VGetArrayPointer(y_sundials);
    realtype* dy_dt = N_VGetArrayPointer(dy_dt_sundials);

#if LOG_ENABLED
    auto& log_state = solver->rhs_stiff_log_state;
    if (!log_state.is_initialized()) {
        log_state.initialize(solver->ySize, solver->process.unitOperations.size());
    }
    log_state.y_copy.assign(y, y + solver->ySize);
#endif

    // Zero-initialize dy_dt before accumulating contributions
    std::fill(dy_dt, dy_dt + solver->ySize, 0.0);

    std::vector<realtype> t_unitOperations;
    t_unitOperations.reserve(solver->process.unitOperations.size());
    for (auto const& unitOperation_ptr : solver->process.unitOperations) {
        size_t first_idx = solver->getUnitOperationStartIdx(unitOperation_ptr.get());
        realtype t_unitOperation = 0.0;
        BENCHMARK(t_unitOperation,
                  { unitOperation_ptr->rhs(t, y + first_idx, dy_dt + first_idx, *solver, true, false, false, false); });
        t_unitOperations.push_back(t_unitOperation);
    }
    // stiff variant has no connections call here; still log components
#if BENCHMARK_ENABLED
    {
        std::ostringstream oss;
        oss.setf(std::ios::scientific);
        oss.precision(6);
        realtype t_total = std::accumulate(t_unitOperations.begin(), t_unitOperations.end(), 0.0);
        oss << "t=" << t << "\t\tt_total=" << (t_total);
        for (size_t i = 0; i < t_unitOperations.size(); ++i) {
            const auto& uo_ptr = solver->process.unitOperations[i];
            oss << "\t\tt_unitOperation[" << i << ",type=" << static_cast<int>(uo_ptr->getType())
                << "]=" << (t_unitOperations[i]);
        }
        oss << "\n";
        LOG_BENCHMARK("rhs_stiff.log", oss.str());
    }
#endif

#if LOG_ENABLED
    log_rhs_call(t, y_sundials, dy_dt_sundials, log_state, "rhs_stiff.log");
#endif

    return 0;
}

#if LOG_ENABLED
// RHS function (static)
int Solver::log_rhs_call(realtype t,
                         N_Vector y_sundials,
                         N_Vector dy_dt_sundials,
                         RhsLogState& state,
                         const std::string& logFileName) {
    // CAREFUL: This function is static, so it cannot access non-static members directly.
    // t and y are not member variables, they are passed as arguments. Access the object variables via solver.t and
    // solver.y user_data points to the Solver instance

    const realtype* y = N_VGetArrayPointer(y_sundials);
    realtype* dy_dt = N_VGetArrayPointer(dy_dt_sundials);
    auto& call_count = state.call_count;
    auto& last_t = state.last_t;
    auto& y_copy = state.y_copy;
    auto& last_y = state.last_y;
    auto& last_dy_dt = state.last_dy_dt;
    const auto ySize = state.y_copy.size();

    // General logging with execution times
    std::ostringstream oss;
    oss.precision(4);
    oss << std::scientific;
    if (call_count < LOG_FIRST_N_CALLS || call_count % LOG_EVERY_N_CALLS == 0) {
        oss << call_count << ". rhs function call at t = " << t << " s\n";
    }

    // Recognized automatic differentiation
    if (t == last_t && call_count < 1000) {
        oss << "Automatic differentiation recognized:\n";
        oss.precision(15);
        for (size_t i = 0; i < ySize; ++i) {
            if (y[i] != last_y[i]) {
                oss << "    y    [" << std::setw(3) << i << "] changed from " << std::setw(25) << last_y[i] << " to "
                    << std::setw(25) << y[i] << " by " << std::setw(25) << (y[i] - last_y[i]) << "\n";
            }
        }
        for (size_t i = 0; i < ySize; ++i) {
            if (dy_dt[i] != last_dy_dt[i]) {
                oss << "    dy_dt[" << std::setw(3) << i << "] changed from " << std::setw(25) << last_dy_dt[i] << " to "
                    << std::setw(25) << dy_dt[i] << " by " << std::setw(25) << (dy_dt[i] - last_dy_dt[i]) << "\n";
            }
        }
        oss.precision(4);
    } else {
        // Store the current state for next comparison
        last_y.assign(y, y + ySize);
        last_dy_dt.assign(dy_dt, dy_dt + ySize);
    }

    // Check for NaN or Inf in y or dy_dt
    bool nan_or_inf = false;
    for (size_t i = 0; i < ySize; ++i) {
        if (std::isnan(dy_dt[i]) || std::isinf(dy_dt[i])) {
            nan_or_inf = true;
        }
        if (std::isnan(y[i]) || std::isinf(y[i])) {
            nan_or_inf = true;
        }
    }
    if (nan_or_inf) {
        for (size_t i = 0; i < ySize; ++i) {
            oss << "\ty[" << i << "] = " << y[i] << ", dy_dt[" << i << "] = " << dy_dt[i] << "\n";
        }
    }
    LOG(logFileName, oss.str());
    if (nan_or_inf) {
        logger::flush_all_logs();
        logger::close_all_logs();
        throw std::runtime_error("y or dy_dt contains NaN or Inf (call_count=" + std::to_string(call_count) +
                                 ", t=" + std::to_string(t) + ")");
    }

    last_t = t;
    call_count++;
    return 0;
}
#endif

#define SV_CELL_I_COMP_J(y, i, j, n_components) N_VGetArrayPointer(y)[i * n_components + j]
#define SBM_CELL_I_COMP_J(J, i, j, di, dj, n_components) SM_ELEMENT_B(J, i* n_components + j, di * n_components + dj)

// Jacobian function: ONLY FOR SIMPLE ONE CELL WATER REACTION EXAMPLE / NOT IMPLEMENTED GENERALLY!
int Solver::jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
    // CAREFUL: This function is static, so it cannot access non-static members directly.

    // Jacobian is defined as J_ij = d(f_i)/d(y_j)
    // SUNBandMatrix uses column-major storage

    throw std::runtime_error(
        "Warning: Jacobian function is only implemented for simple one cell water reaction example! Comment out this "
        "line to proceed anyway.");

    Solver* solver = static_cast<Solver*>(user_data);

    auto n_components = solver->process.unitOperations[0]->n_components();
    auto n_cells = (solver->process.unitOperations[0]->y_size()) / n_components;

    assert(n_components == 27);
    assert(n_cells == 1);

    for (auto i = 0; i < n_cells; i++) {
        auto di = i;  // Reactions only inside cell

        realtype k_f = 1.0;
        realtype k_b = 5.550844e+15;

        auto dr_dH2O = k_f;
        auto dr_dH_plus = -k_b * SV_CELL_I_COMP_J(y, i, 2, n_components);
        auto dr_dOH_minus = -k_b * SV_CELL_I_COMP_J(y, i, 1, n_components);

        SBM_CELL_I_COMP_J(J, i, 0, di, 0, n_components) = -dr_dH2O;
        SBM_CELL_I_COMP_J(J, i, 0, di, 1, n_components) = -dr_dH_plus;
        SBM_CELL_I_COMP_J(J, i, 0, di, 2, n_components) = -dr_dOH_minus;

        SBM_CELL_I_COMP_J(J, i, 1, di, 0, n_components) = dr_dH2O;
        SBM_CELL_I_COMP_J(J, i, 1, di, 1, n_components) = dr_dH_plus;
        SBM_CELL_I_COMP_J(J, i, 1, di, 2, n_components) = dr_dOH_minus;

        SBM_CELL_I_COMP_J(J, i, 1, di, 0, n_components) = dr_dH2O;
        SBM_CELL_I_COMP_J(J, i, 1, di, 1, n_components) = dr_dH_plus;
        SBM_CELL_I_COMP_J(J, i, 1, di, 2, n_components) = dr_dOH_minus;

        LOG("jac.log", std::string("Jacobian call: k_f=") + std::to_string(k_f) + ", k_b=" + std::to_string(k_b) +
                           ", dr_dH2O=" + std::to_string(dr_dH2O) + ", dr_dH_plus=" + std::to_string(dr_dH_plus) +
                           ", dr_dOH_minus=" + std::to_string(dr_dOH_minus) +
                           ", [H+]=" + std::to_string(SV_CELL_I_COMP_J(y, i, 1, n_components)) +
                           ", [OH-]=" + std::to_string(SV_CELL_I_COMP_J(y, i, 2, n_components)) + "\n");
    }

    return 0;
}
