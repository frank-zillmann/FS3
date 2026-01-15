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
#include <stdexcept>
#include <string>

#include "Components/Component.hpp"
#include "EigenDataTypes.hpp"
#include "Logger.hpp"
#include "Process.hpp"
#include "Solver.hpp"

Solver::Solver(const Process& process, SolverType solverType)
    : process(process),
      unitOperationStartIdx(),
      solverType(solverType),
      solver_memory(nullptr),
      y(nullptr),
      J(nullptr),
      lin_sol(nullptr),
      t(0.0),
      ySize(0),
      sunctx(nullptr) {
    // Fill unitOperationStartIdx with the starting index of each unit operation in the solution vector y
    for (const auto& unitOperation_ptr : process.unitOperations) {
        unitOperationStartIdx[unitOperation_ptr.get()] = ySize;
        LOG("solver_initialization.log", "Unit operation " << static_cast<int>(unitOperation_ptr->getType())
                                                           << " starts at index " << ySize << " with y_size "
                                                           << unitOperation_ptr->y_size() << "\n");
        ySize += unitOperation_ptr->y_size();
    }

    // Create Sundials context (must be first)
    int flag = SUNContext_Create(nullptr, &sunctx);
    CHECK_SUNDIALS_FLAG(flag, "SUNContext_Create");

    // Allocate solution vector
    y = N_VNew_Serial(ySize, sunctx);
    if (!y) {
        throw std::runtime_error("Failed to allocate N_Vector y");
    }

    // Set initial conditions from unit operations
    realtype* y_data = N_VGetArrayPointer(y);
    std::fill(y_data, y_data + ySize, 0.0);  // Zero-initialize first

    for (const auto& unitOperation_ptr : process.unitOperations) {
        auto startIdx = getUnitOperationStartIdx(unitOperation_ptr.get());
        if (unitOperation_ptr->y.size() == 0) {
            // Already zero-initialized above -> do nothing
        } else if (unitOperation_ptr->y.size() == unitOperation_ptr->y_size()) {
            ArrayMap mapped_y(y_data + startIdx, unitOperation_ptr->y.rows(), unitOperation_ptr->y.cols(),
                              PhaseStride(unitOperation_ptr->y.cols(), 1));
            mapped_y = unitOperation_ptr->y;
        } else {
            throw std::runtime_error("Unit operation y vector size does not match expected size.");
        }
    }

    // Create solver memory based on selected solver type
    if (solverType == SolverType::BDF) {
        solver_memory = CVodeCreate(CV_BDF, sunctx);
        if (!solver_memory) {
            throw std::runtime_error("Failed to create CVODE BDF solver");
        }
        LOG("solver_initialization.log", "Created CVODE BDF solver.\n");
    } else if (solverType == SolverType::ADAMS) {
        solver_memory = CVodeCreate(CV_ADAMS, sunctx);
        if (!solver_memory) {
            throw std::runtime_error("Failed to create CVODE ADAMS solver");
        }
        LOG("solver_initialization.log", "Created CVODE ADAMS solver.\n");
    } else if (solverType == SolverType::ERK) {
        solver_memory = ERKStepCreate(rhs, t, y, sunctx);
        if (!solver_memory) {
            throw std::runtime_error("Failed to create ERKStep solver");
        }
        LOG("solver_initialization.log", "Created ERKStep solver.\n");
    } else if (solverType == SolverType::ARK) {
        solver_memory = ARKStepCreate(rhs_nonStiff, rhs_stiff, t, y, sunctx);
        if (!solver_memory) {
            throw std::runtime_error("Failed to create ARKStep solver");
        }
        LOG("solver_initialization.log", "Created ARKStep solver.\n");
    } else {
        throw std::invalid_argument("Unsupported SolverType");
    }

    // Further function calls, parameter settings, etc. based on solver type
    if (solverType == SolverType::BDF || solverType == SolverType::ADAMS) {
        // Initialize CVODE memory with RHS function and initial conditions
        flag = CVodeInit(solver_memory, rhs, t, y);
        CHECK_SUNDIALS_FLAG(flag, "CVodeInit");

        // Specify scalar tolerances
        flag = CVodeSStolerances(solver_memory, reltol, abstol);
        CHECK_SUNDIALS_FLAG(flag, "CVodeSStolerances");

        if (use_sundials_non_negative_constraint) {
            constraints = N_VNew_Serial(ySize, sunctx);
            if (!constraints) {
                throw std::runtime_error("Failed to allocate N_Vector constraints");
            }
            realtype* cdata = NV_DATA_S(constraints);
            for (int i = 0; i < ySize; ++i) {
                // 1.0 => enforce y[i] >= 0.0
                // 0.0 => no constraint
                // -1.0 => enforce y[i] <= 0.0
                cdata[i] = 1.0;
            }
            // Register constraints with CVODE (do this after CVodeCreate/CVodeInit, before CVode)
            flag = CVodeSetConstraints(solver_memory, constraints);
            CHECK_SUNDIALS_FLAG(flag, "CVodeSetConstraints");
        }

        if (process.t_end < std::numeric_limits<realtype>::infinity()) {
            // Set the stop time for the integration
            flag = CVodeSetStopTime(solver_memory, process.t_end);
            CHECK_SUNDIALS_FLAG(flag, "CVodeSetStopTime");
        }

        if (solverType == SolverType::BDF) {
            // BDF: uses Newton iteration and therefore requires a linear solver; allow optional Jacobian
            if (use_analytic_jacobian) {
                throw std::runtime_error("Analytic jacobian for BDF is not implemented");
            } else {
                // No user Jacobian: use dense solver; CVODE will approximate Jacobian internally
                J = SUNDenseMatrix(ySize, ySize, sunctx);
                if (!J) {
                    throw std::runtime_error("Failed to create SUNDenseMatrix J");
                }
                lin_sol = SUNLinSol_Dense(y, J, sunctx);
                if (!lin_sol) {
                    throw std::runtime_error("Failed to create SUNLinSol_Dense linear solver");
                }
                flag = CVodeSetLinearSolver(solver_memory, lin_sol, J);
                CHECK_SUNDIALS_FLAG(flag, "CVodeSetLinearSolver");
            }

            // Limit order for stability in stiff bursts (default 5 for BDF)
            flag = CVodeSetMaxOrd(solver_memory, max_order_bdf);
            CHECK_SUNDIALS_FLAG(flag, "CVodeSetMaxOrd");

            // Max nonlinear iterations only relevant for Newton/BDF
            flag = CVodeSetMaxNonlinIters(solver_memory, max_nonlin_solver_iters);
            CHECK_SUNDIALS_FLAG(flag, "CVodeSetMaxNonlinIters");
        } else {
            // ADAMS: use functional (fixed-point) iteration; no linear solver/Jacobian needed
            // Create and attach the fixed-point nonlinear solver
            nls = SUNNonlinSol_FixedPoint(y, 0, sunctx);
            if (!nls) {
                throw std::runtime_error("Failed to create SUNNonlinSol_FixedPoint");
            }
            flag = CVodeSetNonlinearSolver(solver_memory, nls);
            CHECK_SUNDIALS_FLAG(flag, "CVodeSetNonlinearSolver");

            // Ensure we do not hold matrix/linear solver
            J = nullptr;
            lin_sol = nullptr;

            // Limit order for Adams to improve robustness in rapidly varying regions
            flag = CVodeSetMaxOrd(solver_memory, max_order_adams);
            CHECK_SUNDIALS_FLAG(flag, "CVodeSetMaxOrd");
        }

        // Set user data pointer to this instance (for callbacks)
        flag = CVodeSetUserData(solver_memory, this);
        CHECK_SUNDIALS_FLAG(flag, "CVodeSetUserData");

        flag = CVodeSetMaxNumSteps(solver_memory, max_steps);
        CHECK_SUNDIALS_FLAG(flag, "CVodeSetMaxNumSteps");

        // Optional initial step size hint
        if (init_step > 0.0) {
            flag = CVodeSetInitStep(solver_memory, init_step);
            CHECK_SUNDIALS_FLAG(flag, "CVodeSetInitStep");
        }

        // Set maximum allowed step size if user provided a positive cap (SUNDIALS default 0.0 = no cap)
        if (max_step > 0.0) {
            flag = CVodeSetMaxStep(solver_memory, max_step);
            CHECK_SUNDIALS_FLAG(flag, "CVodeSetMaxStep");
        }

        // Allow smaller steps and more recoveries in stiff regions
        flag = CVodeSetMinStep(solver_memory, min_step);
        CHECK_SUNDIALS_FLAG(flag, "CVodeSetMinStep");
        if (adams_stability_limit_detection && solverType == SolverType::ADAMS) {
            flag = CVodeSetStabLimDet(solver_memory, SUNTRUE);
            CHECK_SUNDIALS_FLAG(flag, "CVodeSetStabLimDet");
        }
        flag = CVodeSetMaxErrTestFails(solver_memory, max_err_test_fails);
        CHECK_SUNDIALS_FLAG(flag, "CVodeSetMaxErrTestFails");

    } else if (solverType == SolverType::ERK) {
        // Specify scalar tolerances
        flag = ERKStepSStolerances(solver_memory, reltol, abstol);
        CHECK_SUNDIALS_FLAG(flag, "ERKStepSStolerances");

        flag = ERKStepSetUserData(solver_memory, this);
        CHECK_SUNDIALS_FLAG(flag, "ERKStepSetUserData");

        flag = ERKStepSetMaxNumSteps(solver_memory, max_steps);
        CHECK_SUNDIALS_FLAG(flag, "ERKStepSetMaxNumSteps");

        // // Allow smaller steps if needed near sharp transients/discontinuities
        // flag = ERKStepSetMinStep(solver_memory, min_step);
        // CHECK_SUNDIALS_FLAG(flag, "ERKStepSetMinStep");

        if (process.t_end < std::numeric_limits<realtype>::infinity()) {
            // Set the stop time for the integration
            flag = ERKStepSetStopTime(solver_memory, process.t_end);
            CHECK_SUNDIALS_FLAG(flag, "ERKStepSetStopTime");
        }

        // Robustness settings for explicit solver
        flag = ERKStepSetMinStep(solver_memory, min_step);
        CHECK_SUNDIALS_FLAG(flag, "ERKStepSetMinStep");
        if (max_step > 0.0) {
            flag = ERKStepSetMaxStep(solver_memory, max_step);
            CHECK_SUNDIALS_FLAG(flag, "ERKStepSetMaxStep");
        }
        if (fixed_step > 0.0) {
            // Force constant step size; disables adaptivity
            flag = ERKStepSetFixedStep(solver_memory, fixed_step);
            CHECK_SUNDIALS_FLAG(flag, "ERKStepSetFixedStep");
        }
        flag = ERKStepSetMaxErrTestFails(solver_memory, max_err_test_fails);
        CHECK_SUNDIALS_FLAG(flag, "ERKStepSetMaxErrTestFails");

    } else if (solverType == SolverType::ARK) {
        // Specify scalar tolerances
        flag = ARKStepSStolerances(solver_memory, reltol, abstol);
        CHECK_SUNDIALS_FLAG(flag, "ARKStepSStolerances");

        flag = ARKStepSetUserData(solver_memory, this);
        CHECK_SUNDIALS_FLAG(flag, "ARKStepSetUserData");

        flag = ARKStepSetMaxNumSteps(solver_memory, max_steps);
        CHECK_SUNDIALS_FLAG(flag, "ARKStepSetMaxNumSteps");

        if (process.t_end < std::numeric_limits<realtype>::infinity()) {
            // Set the stop time for the integration
            flag = ARKStepSetStopTime(solver_memory, process.t_end);
            CHECK_SUNDIALS_FLAG(flag, "ARKStepSetStopTime");
        }

        // ARKode with implicit/stiff part always needs a linear solver
        auto half_band_width = process.componentSystem.n_components;
        J = SUNBandMatrix(ySize, half_band_width, half_band_width, sunctx);
        if (!J) {
            throw std::runtime_error("Failed to create SUNBandMatrix J");
        }
        lin_sol = SUNLinSol_Band(y, J, sunctx);
        if (!lin_sol) {
            throw std::runtime_error("Failed to create SUNLinSol_Band linear solver");
        }

        flag = ARKStepSetLinearSolver(solver_memory, lin_sol, J);
        CHECK_SUNDIALS_FLAG(flag, "ARKStepSetLinearSolver");

        if (use_analytic_jacobian) {
            std::cout << "Jacobian is used!";
            // Set the Jacobian function callback for Newton iterations
            flag = ARKStepSetJacFn(solver_memory, jac);
            CHECK_SUNDIALS_FLAG(flag, "ARKStepSetJacFn");
        }
        // Note: If use_jacobian is false, ARKode will use finite difference approximation

        // Robustness: permit very small steps and more error test recoveries
        flag = ARKStepSetMinStep(solver_memory, min_step);
        CHECK_SUNDIALS_FLAG(flag, "ARKStepSetMinStep");
        if (max_step > 0.0) {
            flag = ARKStepSetMaxStep(solver_memory, max_step);
            CHECK_SUNDIALS_FLAG(flag, "ARKStepSetMaxStep");
        }
        if (fixed_step > 0.0) {
            flag = ARKStepSetFixedStep(solver_memory, fixed_step);
            CHECK_SUNDIALS_FLAG(flag, "ARKStepSetFixedStep");
        }
        flag = ARKStepSetMaxErrTestFails(solver_memory, max_err_test_fails);
        CHECK_SUNDIALS_FLAG(flag, "ARKStepSetMaxErrTestFails");
    } else {
        throw std::invalid_argument("Unsupported SolverType");
    }
}

// Destructor cleans up all allocated Sundials objects
Solver::~Solver() {
    if (solver_memory && (solverType == SolverType::ADAMS || solverType == SolverType::BDF)) CVodeFree(&solver_memory);
    if (solver_memory && solverType == SolverType::ERK) ERKStepFree(&solver_memory);
    if (solver_memory && solverType == SolverType::ARK) ARKStepFree(&solver_memory);

    if (y) N_VDestroy(y);
    if (constraints) N_VDestroy(constraints);
    if (J) SUNMatDestroy(J);
    if (lin_sol) SUNLinSolFree(lin_sol);
    if (nls) SUNNonlinSolFree(nls);
    if (sunctx) SUNContext_Free(&sunctx);
}

const std::vector<realtype> Solver::getY() const {
    std::vector<realtype> y_vec;
    if (y != nullptr && ySize > 0) {
        y_vec.reserve(ySize);
        const realtype* y_data = N_VGetArrayPointer(y);
        y_vec.assign(y_data, y_data + ySize);
    }
    return y_vec;
}

std::size_t Solver::getUnitOperationStartIdx(const UnitOperationBase* unitOp) const {
    auto it = unitOperationStartIdx.find(unitOp);
    if (it != unitOperationStartIdx.end()) {
        return it->second;
    }
    throw std::out_of_range("UnitOperationBase not found in unitOperationStartIdx map.");
}
