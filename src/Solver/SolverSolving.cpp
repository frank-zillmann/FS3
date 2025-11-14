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
#include "cnpy.h"

bool Solver::checkTimeout() const {
    // Check timeout periodically (not every step to avoid overhead)
    auto elapsed = std::chrono::steady_clock::now() - t_solve_start;
    elapsed_seconds = std::chrono::duration<realtype>(elapsed).count();

    if (elapsed_seconds >= timeout_seconds) {
        return true;
    } else {
        return false;
    }
}

int Solver::run_solver(realtype t_stop) {
    int flag = 0;
    if (internal_time_stamps.empty()) internal_time_stamps.push_back(t);

    if (solverType == SolverType::BDF || solverType == SolverType::ADAMS) {
        CVodeSetStopTime(solver_memory, t_stop);
    } else if (solverType == SolverType::ERK) {
        ERKStepSetStopTime(solver_memory, t_stop);
    } else if (solverType == SolverType::ARK) {
        ARKStepSetStopTime(solver_memory, t_stop);
    } else {
        throw std::invalid_argument("Unsupported SolverType");
    }

    while (t < t_stop) {
        realtype t_step = 0.0;  // BENCHMARK measures seconds
        BENCHMARK(t_step, {
            if (solverType == SolverType::BDF || solverType == SolverType::ADAMS) {
                flag = CVode(solver_memory, t_stop, y, &t, CV_ONE_STEP);
            } else if (solverType == SolverType::ERK) {
                flag = ERKStepEvolve(solver_memory, t_stop, y, &t, ARK_ONE_STEP);
            } else if (solverType == SolverType::ARK) {
                flag = ARKStepEvolve(solver_memory, t_stop, y, &t, ARK_ONE_STEP);
            } else {
                throw std::invalid_argument("Unsupported SolverType");
            }
        });
#if BENCHMARK_ENABLED
        {
            std::ostringstream oss;
            oss.setf(std::ios::scientific);
            oss.precision(6);
            oss << "step=" << (n_steps + 1) << "\t\tt=" << t << "\t\tt_step=" << (t_step) << "\n";
            LOG_BENCHMARK("solver_step.log", oss.str());
        }
#endif

        n_steps++;
        internal_time_stamps.push_back(t);

        if (flag < 0) {
            logStatistics();
            logger::flush_all_logs();
            logger::close_all_logs();
            CHECK_SUNDIALS_FLAG(flag, "Evolve");
        }
        if (n_steps >= max_steps) {
            logStatistics();
            logger::flush_all_logs();
            logger::close_all_logs();
            std::cerr << "Maximum number of steps (" << max_steps << ") exceeded." << std::endl;
            throw std::runtime_error("Maximum number of steps (" + std::to_string(max_steps) + ") exceeded.");
        }
        steps_since_timeout_check++;
        if (steps_since_timeout_check >= timeout_check_interval) {
            if (checkTimeout()) {
                logStatistics();
                logger::flush_all_logs();
                logger::close_all_logs();
                std::cerr << "Solver timed out after " << elapsed_seconds << " seconds." << std::endl;
                throw std::runtime_error("Solver timed out after " + std::to_string(elapsed_seconds) + " seconds.");
            }
            steps_since_timeout_check = 0;
        }
    }
    return flag;
}

// Simple solve step: integrate to time tout
void Solver::solve(realtype t_stop, realtype maxSolveTime) {
    // Initialize timeout tracking
    this->timeout_seconds = maxSolveTime;
    this->t_solve_start = std::chrono::steady_clock::now();
    this->steps_since_timeout_check = 0;

    int flag = 0;
    realtype* y_ptr = N_VGetArrayPointer(y);

    // Sort Observers by time before processing
    std::sort(observers.begin(), observers.end(),
              [](const std::shared_ptr<SnapshotObserver>& a, const std::shared_ptr<SnapshotObserver>& b) {
                  return a->t_desired < b->t_desired;
              });

    // Iterate over snapshot observers and write snapshots at specified times
    for (auto& observer : observers) {
        const auto snapshotTime = observer->t_desired;
        if (snapshotTime > t_stop) break;

        if (snapshotTime > t) {
            flag = run_solver(snapshotTime);
        }

        realtype max_time_diff = 0.1;
        if (std::abs(t - snapshotTime) > max_time_diff) {
            throw std::runtime_error("Current time t=" + std::to_string(t) + " differs by more than " +
                                     std::to_string(max_time_diff) + " from snapshot time " +
                                     std::to_string(snapshotTime) + "!");
        }

        observer->write(t, y_ptr + getUnitOperationStartIdx(observer->mapper.unitOperation_ptr));
    }

    flag = run_solver(t_stop);

    checkTimeout();
    logStatistics();
}

void Solver::logStatistics() {
#if LOG_ENABLED
    if (solverType == SolverType::BDF || solverType == SolverType::ADAMS) {
        std::ostringstream oss;
        oss << std::setprecision(6);

        // Basic CVODE counters
        long int nsteps = 0, nfevals = 0, netfails = 0, nlinsetups = 0;
        if (CVodeGetNumSteps(solver_memory, &nsteps) == CV_SUCCESS) {
            oss << "nsteps=" << nsteps << "  ";
        }
        if (CVodeGetNumRhsEvals(solver_memory, &nfevals) == CV_SUCCESS) {
            oss << "nfevals=" << nfevals << "  ";
        }
        if (CVodeGetNumErrTestFails(solver_memory, &netfails) == CV_SUCCESS) {
            oss << "errtestfails=" << netfails << "  ";
        }
        if (CVodeGetNumLinSolvSetups(solver_memory, &nlinsetups) == CV_SUCCESS) {
            oss << "linSetups=" << nlinsetups << "  ";
        }

        // Nonlinear solver statistics
        long int nniters = 0, nnconvfails = 0;
        if (CVodeGetNumNonlinSolvIters(solver_memory, &nniters) == CV_SUCCESS) {
            oss << "nonlinIters=" << nniters << "  ";
        }
        if (CVodeGetNumNonlinSolvConvFails(solver_memory, &nnconvfails) == CV_SUCCESS) {
            oss << "nonlinConvFails=" << nnconvfails << "  ";
        }

        // Jacobian / Jtimes counts only when a linear solver is present (e.g., BDF)
        if (solverType == SolverType::BDF) {
            long int njacevals = 0, njtimes = 0;
            if (CVodeGetNumJacEvals(solver_memory, &njacevals) == CV_SUCCESS) {
                oss << "JacEvals=" << njacevals << "  ";
            }
            if (CVodeGetNumJtimesEvals(solver_memory, &njtimes) == CV_SUCCESS) {
                oss << "JtimesEvals=" << njtimes << "  ";
            }
        }

        // Step / order / t / gamma diagnostics
        realtype last_h = 0.0, cur_h = 0.0, tcur = 0.0, gamma = 0.0;
        int qcur = 0, qlast = 0;
        if (CVodeGetLastStep(solver_memory, &last_h) == CV_SUCCESS)
            oss << "last_h=" << std::scientific << last_h << "  ";
        if (CVodeGetCurrentStep(solver_memory, &cur_h) == CV_SUCCESS)
            oss << "cur_h=" << std::scientific << cur_h << "  ";
        if (CVodeGetCurrentTime(solver_memory, &tcur) == CV_SUCCESS) oss << "tcur=" << tcur << "  ";
        if (CVodeGetCurrentOrder(solver_memory, &qcur) == CV_SUCCESS) oss << "order=" << qcur << "  ";
        if (CVodeGetLastOrder(solver_memory, &qlast) == CV_SUCCESS) oss << "last_order=" << qlast << "  ";
        if (CVodeGetCurrentGamma(solver_memory, &gamma) == CV_SUCCESS)
            oss << "gamma=" << std::scientific << gamma << "  ";

        oss << "\n";

        LOG("solver_statistics.log", oss.str());

    } else if (solverType == SolverType::ERK) {
        // long int nsteps = 0, nfevals = 0, netfails = 0;
        // if (ERKStepGetNumSteps(solver_memory, &nsteps) == CV_SUCCESS) {
        //     oss << "nsteps=" << nsteps << "  ";
        // }
        // if (ERKStepGetNumRhsEvals(solver_memory, &nfevals) == CV_SUCCESS) {
        //     oss << "nfevals=" << nfevals << "  ";
        // }
        // if (ERKStepGetNumErrTestFails(solver_memory, &netfails) == CV_SUCCESS) {
        //     oss << "errtestfails=" << netfails << "  ";
        // }

        // // Step / order / t diagnostics
        // realtype last_h = 0.0, cur_h = 0.0, tcur = 0.0;
        // int qcur = 0;
        // if (ERKStepGetLastStep(solver_memory, &last_h) == CV_SUCCESS)
        //     oss << "last_h=" << std::scientific << last_h << "  ";
        // if (ERKStepGetCurrentStep(solver_memory, &cur_h) == CV_SUCCESS)
        //     oss << "cur_h=" << std::scientific << cur_h << "  ";
        // if (ERKStepGetCurrentTime(solver_memory, &tcur) == CV_SUCCESS) oss << "tcur=" << tcur << "  ";
        // if (ERKStepGetCurrentOrder(solver_memory, &qcur) == CV_SUCCESS) oss << "order=" << qcur << "  ";

        std::string path = logger::log_dir() + "/solver_statistics.log";
        FILE* fp = fopen(path.c_str(), "a");
        if (!fp) {
            LOG("solver_statistics.log", "Failed to open log file for ERKStepWriteParameters");
            return;
        }

        ERKStepWriteButcher(solver_memory, fp);
        ERKStepWriteParameters(solver_memory, fp);

        fclose(fp);  // safe since it’s not the logger’s stream

    } else if (solverType == SolverType::ARK) {
        std::ostringstream oss;
        oss << std::setprecision(6);

        // Basic ARKStep counters
        long int nsteps = 0, nfevals_explicit = 0, nfevals_implicit = 0, netfails = 0, nlinsetups = 0;
        if (ARKStepGetNumSteps(solver_memory, &nsteps) == ARK_SUCCESS) {
            oss << "nsteps=" << nsteps << "  ";
        }
        if (ARKStepGetNumRhsEvals(solver_memory, &nfevals_explicit, &nfevals_implicit) == ARK_SUCCESS) {
            oss << "nfevals_explicit=" << nfevals_explicit << "  ";
            oss << "nfevals_implicit=" << nfevals_implicit << "  ";
        }
        if (ARKStepGetNumErrTestFails(solver_memory, &netfails) == ARK_SUCCESS) {
            oss << "errtestfails=" << netfails << "  ";
        }
        if (ARKStepGetNumLinSolvSetups(solver_memory, &nlinsetups) == ARK_SUCCESS) {
            oss << "linSetups=" << nlinsetups << "  ";
        }

        // Nonlinear solver statistics
        long int nniters = 0, nnconvfails = 0;
        if (ARKStepGetNumNonlinSolvIters(solver_memory, &nniters) == ARK_SUCCESS) {
            oss << "nonlinIters=" << nniters << "  ";
        }
        if (ARKStepGetNumNonlinSolvConvFails(solver_memory, &nnconvfails) == ARK_SUCCESS) {
            oss << "nonlinConvFails=" << nnconvfails << "  ";
        }

        // Jacobian evaluations (if available)
        long int njacevals = 0;
        if (ARKStepGetNumJacEvals(solver_memory, &njacevals) == ARK_SUCCESS) {
            oss << "JacEvals=" << njacevals << "  ";
        }

        // Step / t diagnostics
        realtype last_h = 0.0, cur_h = 0.0, tcur = 0.0;
        if (ARKStepGetLastStep(solver_memory, &last_h) == ARK_SUCCESS)
            oss << "last_h=" << std::scientific << last_h << "  ";
        if (ARKStepGetCurrentStep(solver_memory, &cur_h) == ARK_SUCCESS)
            oss << "cur_h=" << std::scientific << cur_h << "  ";
        if (ARKStepGetCurrentTime(solver_memory, &tcur) == ARK_SUCCESS) oss << "tcur=" << tcur << "  ";

        oss << "\n";

        LOG("solver_statistics.log", oss.str());

    } else {
        throw std::invalid_argument("Unsupported SolverType");
    }
#endif
}
