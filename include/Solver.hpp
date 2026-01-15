#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_erkstep.h>
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>
// Nonlinear solvers
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>

#include <chrono>
#include <ctime>
#include <limits>
#include <memory>
#include <unordered_map>
#include <vector>

#include "EigenDataTypes.hpp"
#include "Logger.hpp"
#include "Observers/SnapshotObserver.hpp"

// Forward declarations
class Process;
class UnitOperationBase;

/// @brief Available solver types for ODE integration
enum class SolverType { BDF, ADAMS, ERK, ARK };

/**
 * @brief SUNDIALS-based ODE solver wrapper
 *
 * Provides a high-level interface to CVODE/ARKODE (BDF, Adams, ERK, ARK),
 * builds global indexing, configures tolerances/limits, and exposes callbacks
 * for RHS and optional Jacobian. Logs internal time stamps and statistics.
 */
class Solver {
   public:
    Solver(const Process& process, SolverType solverType = SolverType::BDF);
    ~Solver();

    void solve(realtype t_stop = std::numeric_limits<realtype>::infinity(),
               realtype timeout_seconds = std::numeric_limits<realtype>::infinity());

    int run_solver(realtype t_stop);

    sunindextype getYSize() const { return ySize; }
    const std::vector<realtype> getY() const;
    const realtype getT() const { return t; }
    const std::vector<realtype>& getInternalTimeStamps() const { return internal_time_stamps; }

    std::size_t getUnitOperationStartIdx(const UnitOperationBase* unitOp) const;
    void logStatistics();

    void add(std::shared_ptr<SnapshotObserver> observer) { observers.push_back(observer); }

    // Get the elapsed wall-clock time of the last solve() call in seconds
    realtype getSolveTime() const { return elapsed_seconds; }

   protected:
    // RHS function f(t,y) = y'
    static int rhs(realtype t, N_Vector y, N_Vector dy_dt, void* user_data);
    static int rhs_stiff(realtype t, N_Vector y_sundials, N_Vector dy_dt_sundials, void* user_data);
    static int rhs_nonStiff(realtype t, N_Vector y_sundials, N_Vector dy_dt_sundials, void* user_data);

    // Jacobian function
    static int jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

    bool checkTimeout() const;

    const Process& process;
    std::unordered_map<const UnitOperationBase*, std::size_t> unitOperationStartIdx;

    // Type of the solver
    SolverType solverType;

    // Simulation parameters
    realtype reltol = 1e-5;
    realtype abstol = 1e-12;

    // Step size and solver robustness controls
    // h_min: minimum allowed step size. SUNDIALS default: 0.0 (no lower bound) for CVODE/ERK/ARK
    realtype min_step = 0.0;
    // h_max: maximum allowed step size. SUNDIALS default: 0.0 (no upper bound) for CVODE/ERK/ARK
    realtype max_step = 0.0;
    // Error test failures allowed per (macro-)step. SUNDIALS default: 7 for CVODE/ERK/ARK
    long int max_err_test_fails = 20;  // allow more recoveries through step rejections than default
    // Max method order. SUNDIALS defaults: BDF=5, Adams=12
    int max_order_bdf = 5;     // CVODE BDF default: 5
    int max_order_adams = 12;  // CVODE Adams default: 12
    // Max nonlinear iterations per step (Newton/fixed-point). SUNDIALS CVODE default: 3
    int max_nonlin_solver_iters = 10;
    // Use user-supplied analytic Jacobian (otherwise finite-difference/JTimes). Default in SUNDIALS: false
    bool use_analytic_jacobian = false;
    // Initial step size hint. SUNDIALS default: 0.0 (internally estimated).
    realtype init_step = 0.0;
    // Adams stability limit detection. SUNDIALS default: false
    // TODO: Try to turn this on
    bool adams_stability_limit_detection = false;
    // Fixed step size (only for ERK/ARK). SUNDIALS default: disabled; 0.0 means adaptive.
    realtype fixed_step = 0.0;

    // Constraints (tells Sundials y >= 0 for all times)
    bool use_sundials_non_negative_constraint = false;
    N_Vector constraints = nullptr;

    // Sundials objects
    SUNContext sunctx;
    void* solver_memory;
    SUNLinearSolver lin_sol;
    SUNNonlinearSolver nls = nullptr;

    // State variables
    N_Vector y;
    SUNMatrix J;
    realtype t;
    sunindextype ySize;

    // Observers
    std::vector<std::shared_ptr<SnapshotObserver>> observers;

    long int max_steps = std::numeric_limits<long int>::max();
    long int n_steps = 0;

    // Stores all internal time stamps chosen by Sundials during integration
    std::vector<realtype> internal_time_stamps;

    // Timeout functionality
    mutable realtype timeout_seconds = std::numeric_limits<realtype>::infinity();  // No timeout by default
    mutable std::chrono::steady_clock::time_point t_solve_start;
    mutable int timeout_check_interval = 100;  // Check timeout every N steps to avoid overhead
    mutable int steps_since_timeout_check = 0;
    mutable realtype elapsed_seconds;  // Time elapsed since start of solve() in seconds

#if LOG_ENABLED
    // Helper struct to track state for each RHS method independently
    struct RhsLogState {
        int call_count = 1;
        realtype last_t = -1.0;
        std::vector<realtype> y_copy;
        std::vector<realtype> last_y;
        std::vector<realtype> last_dy_dt;

        RhsLogState() = default;

        void initialize(size_t y_size, size_t n_unit_operations) {
            y_copy.resize(y_size, 0.0);
            last_y.resize(y_size, 0.0);
            last_dy_dt.resize(y_size, 0.0);
        }

        bool is_initialized() const {
            if (y_copy.empty() || last_y.empty() || last_dy_dt.empty()) {
                return false;
            }
            return true;
        }
    };

    // Common implementation for all RHS methods
    static int log_rhs_call(realtype t,
                            N_Vector y_sundials,
                            N_Vector dy_dt_sundials,
                            RhsLogState& state,
                            const std::string& logFileName);

    // Separate tracking state for each RHS method
    RhsLogState rhs_log_state;
    RhsLogState rhs_stiff_log_state;
    RhsLogState rhs_nonStiff_log_state;
#endif
};

// Helper macro for error checking
#define CHECK_SUNDIALS_FLAG(flag, funcname)                                                              \
    if ((flag) != CV_SUCCESS) {                                                                          \
        throw std::runtime_error(std::string(funcname) + " failed with flag = " + std::to_string(flag)); \
    }

#endif  // SOLVER_HPP