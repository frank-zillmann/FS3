#include "UnitOperations/RS_MagneticCaptureProcessChamber.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <string>

#include "ConvectionDispersionOperators.hpp"
#include "EigenDataTypes.hpp"
#include "Logger.hpp"
#include "PhysicalConstants.hpp"

sunindextype RS_MagneticCaptureProcessChamber::idx(sunindextype finiteVolumeIdx, sunindextype componentIdx, Phase phase) const {
    return UnitOperationBase::idx<Phase::Liquid, Phase::SolidOrSlurry>(finiteVolumeIdx, componentIdx, phase);
}

void RS_MagneticCaptureProcessChamber::init() {
    // Find ranges of different types of components
    const auto& componentSystem = reactionSystem.componentSystem;
    const auto& components = componentSystem.components;

    NMC_range = {0, componentSystem.getIdx("hIgG")};
    MNPG_range = {componentSystem.getIdx("MNP-OH₂⁺"), componentSystem.getIdx("MNP-hIgG")};
    MNP_range = {componentSystem.getIdx("MNP1"), componentSystem.getIdx("MNP10")};

    // Resize arrays based on number of magnetic components
    auto n_MNP_components = MNP_range.second - MNP_range.first + 1;
    radii_MNP.resize(n_MNP_components);
    densities_MNP.resize(n_MNP_components);
    magnetic_saturations_MNP.resize(n_MNP_components);

    // Second pass: fill the arrays with values
    sunindextype vectorIndex = 0;
    for (sunindextype componentIndex = MNP_range.first; componentIndex <= MNP_range.second; componentIndex++) {
        if (componentSystem.components[componentIndex].radius.has_value()) {
            radii_MNP(0, vectorIndex) = componentSystem.components[componentIndex].radius.value();
        } else {
            throw std::runtime_error("Magnetic component must have a radius defined!");
        }
        if (componentSystem.components[componentIndex].density.has_value()) {
            densities_MNP(0, vectorIndex) = componentSystem.components[componentIndex].density.value();
        } else {
            throw std::runtime_error("Magnetic component must have a density defined!");
        }
        if (componentSystem.components[componentIndex].magnetic_saturation.has_value()) {
            magnetic_saturations_MNP(0,
                                     vectorIndex) = componentSystem.components[componentIndex].magnetic_saturation.value();
        } else {
            throw std::runtime_error("Magnetic component must have a magnetic saturation defined!");
        }
        vectorIndex++;
    }
}

const ColVector& RS_MagneticCaptureProcessChamber::get_V_l(realtype t, const realtype* y) const {
    if ((t == -1.0) && (y == nullptr)) {
        return V_l;  // arguments indicates that last_V_l can be returned
    } else {
        PhaseStride two_phase_stride(2 * n_components(), 1);

        ConstArrayMap m_sl(y + n_components(), n_cells_per_phase, n_components(), two_phase_stride);

        auto m_sl_MNP = m_sl.middleCols(MNP_range.first, MNP_range.second - MNP_range.first + 1);

        // Calculate cell volumes of fluid and slurry
        auto V_l_and_sl = cell_volume_l_and_sl * ColVector::Ones(n_cells_per_phase);  // (n_cells, 1)
        auto V_sl = alpha_fluid_particle_volume_ratio *
                    (m_sl_MNP.rowwise() / densities_MNP).rowwise().sum();  // (n_cells, 1)
        V_l = V_l_and_sl - V_sl;  // (n_cells, 1) // safe into member variable for later use in connections

        return V_l;
    }
};

realtype RS_MagneticCaptureProcessChamber::errorFunction(realtype t, const realtype* y) const {
    PhaseStride two_phase_stride(2 * n_components(), 1);

    ConstArrayMap m_l(y, n_cells_per_phase, n_components(), two_phase_stride);
    ConstArrayMap m_sl(y + n_components(), n_cells_per_phase, n_components(), two_phase_stride);

    auto m_sl_MNP = m_sl.middleCols(MNP_range.first, MNP_range.second - MNP_range.first + 1);

    // Calculate cell volumes
    auto V_l_and_sl = cell_volume_l_and_sl * ColVector::Ones(n_cells_per_phase);
    auto V_sl = alpha_fluid_particle_volume_ratio * (m_sl_MNP.rowwise() / densities_MNP).rowwise().sum();
    V_l = V_l_and_sl - V_sl;

    // Calculate concentrations
    c_l = m_l.colwise() / V_l;
    c_l = c_l.cwiseMax(0);  // Clamp negatives to 0

    c_sl = m_sl.colwise() / V_sl;
    c_sl = c_sl.cwiseMax(0);  // Clamp negatives to 0
    c_sl = (V_sl / cell_volume_l_and_sl < 0.01)
               .rowwise()
               .replicate(n_components())
               .select(Array::Zero(n_cells_per_phase, n_components()), c_sl);

    activities_l.setZero();
    activities_sl.setZero();

    realtype error_l = reactionSystem.errorFunction(t, c_l, activities_l);
    realtype error_sl = reactionSystem.errorFunction(t, c_sl, activities_sl);

    return std::max(error_l, error_sl);
}

void RS_MagneticCaptureProcessChamber::rhs(realtype t,
                                           const realtype* y,
                                           realtype* dy_dt,
                                           const Solver& solver,
                                           bool enable_reactions,
                                           bool enable_convection,
                                           bool enable_dispersion,
                                           bool enable_otherPhysics) {
    PhaseStride two_phase_stride(2 * n_components(), 1);

    ConstArrayMap m_l(y, n_cells_per_phase, n_components(), two_phase_stride);
    ArrayMap dm_dt_l(dy_dt, n_cells_per_phase, n_components(), two_phase_stride);

    ConstArrayMap m_sl(y + n_components(), n_cells_per_phase, n_components(), two_phase_stride);
    ArrayMap dm_dt_sl(dy_dt + n_components(), n_cells_per_phase, n_components(), two_phase_stride);

    // Create convenient views for subsets of components
    auto m_l_NMC = m_l.middleCols(NMC_range.first, NMC_range.second - NMC_range.first + 1);
    auto m_sl_NMC = m_sl.middleCols(NMC_range.first, NMC_range.second - NMC_range.first + 1);
    auto dm_dt_l_NMC = dm_dt_l.middleCols(NMC_range.first, NMC_range.second - NMC_range.first + 1);
    auto dm_dt_sl_NMC = dm_dt_sl.middleCols(NMC_range.first, NMC_range.second - NMC_range.first + 1);

    auto m_l_MNPG = m_l.middleCols(MNPG_range.first, MNPG_range.second - MNPG_range.first + 1);
    auto m_sl_MNPG = m_sl.middleCols(MNPG_range.first, MNPG_range.second - MNPG_range.first + 1);
    auto dm_dt_l_MNPG = dm_dt_l.middleCols(MNPG_range.first, MNPG_range.second - MNPG_range.first + 1);
    auto dm_dt_sl_MNPG = dm_dt_sl.middleCols(MNPG_range.first, MNPG_range.second - MNPG_range.first + 1);

    auto m_l_MNP = m_l.middleCols(MNP_range.first, MNP_range.second - MNP_range.first + 1);
    auto m_sl_MNP = m_sl.middleCols(MNP_range.first, MNP_range.second - MNP_range.first + 1);
    auto dm_dt_l_MNP = dm_dt_l.middleCols(MNP_range.first, MNP_range.second - MNP_range.first + 1);
    auto dm_dt_sl_MNP = dm_dt_sl.middleCols(MNP_range.first, MNP_range.second - MNP_range.first + 1);

    // Calculate cell volumes of fluid and slurry
    auto V_l_and_sl = cell_volume_l_and_sl * ColVector::Ones(n_cells_per_phase);  // (n_cells, 1)
    auto V_sl = alpha_fluid_particle_volume_ratio * (m_sl_MNP.rowwise() / densities_MNP).rowwise().sum();  // (n_cells, 1)
    V_l = V_l_and_sl - V_sl;  // (n_cells, 1) // safe into member variable for later use in connections

    const auto flowRate = flowRateFunction(t);
    assert(flowRate >= 0.0 && "Flow rate must be non-negative!");

    const auto dispersion_coefficient = dispersion_coefficient_function(t);
    assert(dispersion_coefficient >= 0.0 && "Dispersion coefficient must be non-negative!");

    // Capture physics (enable_otherPhysics controls capture/uncapture)
    // TODO: add Washburn model
    realtype t_capture = 0.0;
    BENCHMARK(t_capture, {
        if (enable_otherPhysics) {
            auto capture_value = capture_function(t);

            if (capture_value == 1.0) {
                // u_0 is pseudo velocity in completely empty chamber
                auto u_0 = flowRate / crossSectionArea;

                // saturation magnzetization is fine (no need to change to Magnetization Function in alche paper)
                auto u_m = (2 * constants::vacuum_permeability * radii_MNP * radii_MNP * magnetic_saturations_MNP *
                            magnetic_field_strength) /
                           (9 * reactionSystem.componentSystem.dynamic_viscosity * disk_height);  // (1, n_MNP_components)

                // MNP capture cells
                auto m_sl_MNP_sum = m_sl_MNP.rowwise().sum();  // (n_cells, 1)
                auto m_sl_MNP_max = MNP_capacity / n_cells_per_phase;
                auto G = ColVector::Ones(n_cells_per_phase) -
                         (m_sl_MNP_sum / m_sl_MNP_max).pow(deposition_rate);  // (n_cells, 1)

                auto a_eff = (u_m / u_0).unaryExpr(a_eff_function);  // (1, n_MNP_components)

                // dm_dt_capture_MNP_wrong = (((n_disks * u_0) / length) * ((m_l_MNP.rowwise() * a_eff).colwise() * G)).eval(); // (n_cells, n_MNP_components)

                dm_dt_capture_MNP = ((flowRate * n_disks / static_cast<realtype>(n_cells_per_phase)) *
                                     (((m_l_MNP.rowwise() * a_eff).colwise() * G).colwise() * V_l.cwiseInverse()))
                                        .eval();  // (n_cells, n_MNP_components)

#if LOG_ENABLED
                if (call_count < LOG_FIRST_N_CALLS || call_count % LOG_EVERY_N_CALLS == 0) {
                    auto porosity = V_l.array() / cell_volume_l_and_sl;

                    std::ostringstream oss;
                    std::streamsize oss_default_precision = oss.precision();
                    oss << call_count << ".\tt = " << std::setprecision(std::numeric_limits<double>::max_digits10) << t;
                    oss << std::setprecision(oss_default_precision) << ", flowRate = " << flowRate
                        << ", magnetic_field_strength = " << magnetic_field_strength << "\n";
                    oss << ", porosity = " << porosity.transpose() << "\n";
                    oss << "V_l = " << V_l.transpose() << "\n";
                    oss << "V_sl = " << V_sl.transpose() << "\n";
                    oss << "porosity = " << porosity.transpose() << "\n";
                    oss << "u_m/u_0 = " << u_m / u_0 << "\n";
                    oss << "a_eff = " << a_eff << "\n";
                    oss << "G = " << G.transpose() << "\n";
                    oss << "dm_dt_capture_MNP = \n" << dm_dt_capture_MNP << "\n";
                    LOG("rhs_capture.log", oss.str());
                }
#endif

            } else if (capture_value <= 0.0) {
                dm_dt_capture_MNP = (capture_value * m_sl_MNP).eval();

#if LOG_ENABLED
                if (call_count < LOG_FIRST_N_CALLS || call_count % LOG_EVERY_N_CALLS == 0) {
                    auto porosity = V_l.array() / cell_volume_l_and_sl;

                    std::ostringstream oss;
                    std::streamsize oss_default_precision = oss.precision();
                    oss << call_count << ".\tt = " << std::setprecision(std::numeric_limits<double>::max_digits10) << t;
                    oss << std::setprecision(oss_default_precision);
                    oss << "porosity = " << porosity.transpose() << "\n";
                    oss << "dm_dt_capture_MNP (uncapture)= \n" << dm_dt_capture_MNP << "\n";
                    LOG("rhs_capture.log", oss.str());
                }
#endif
            } else {
                std::cerr << "Warning: capture function returned unexpected value: " << capture_value << " at t = " << t
                          << " ! No capture/uncapture happens\n";
                dm_dt_capture_MNP = Array::Zero(n_cells_per_phase, MNP_range.second - MNP_range.first + 1);
            }

            // Volume change rates
            auto dm_dt_capture_MNP_sum_over_density = (dm_dt_capture_MNP.rowwise() / densities_MNP).rowwise().sum();  // (n_cells, 1)
            auto dV_dt_sl = alpha_fluid_particle_volume_ratio * dm_dt_capture_MNP_sum_over_density;  // (n_cells, 1)
            auto dV_dt_l = -dV_dt_sl;                                                                // (n_cells, 1)

            // NMC capture (Mechanic capture on MNPs and MNPGs neglected)

            // Old only capture case version:
            // auto dm_dt_capture_NMC = (m_l_NMC.colwise() * dV_dt_sl).colwise() / V_l;  // (n_cells, n_NMC_components)

            // For cells where dV_dt_sl < 0, use m_sl_NMC and V_sl instead of m_l_NMC and V_l
            auto mask_capture_uncapture = (dV_dt_sl.array() >= 0).cast<realtype>();  // (n_cells, 1), 1.0 where dV_dt_sl < 0, else 0.0

            auto m_from_NMC = m_l_NMC.colwise() * mask_capture_uncapture +
                              m_sl_NMC.colwise() * (1 - mask_capture_uncapture);
            auto V_from_NMC = V_l.colwise() * mask_capture_uncapture + V_sl.colwise() * (1 - mask_capture_uncapture);

            auto dm_dt_capture_NMC = (m_from_NMC.colwise() * dV_dt_sl).colwise() /
                                     (V_from_NMC + 1e-20);  // (n_cells, n_NMC_components)

            // MNPG capture

            // Old only capture case version:
            // auto m_l_MNP_sum = m_l_MNP.rowwise().sum();                                  // (n_cells, 1)
            // auto MNPG_MNP_ratio = m_l_MNPG.colwise() / (m_l_MNP_sum + 1e-12);            // (n_cells, n_MNPG_components)
            // auto dm_dt_capture_MNP_sum = dm_dt_capture_MNP.rowwise().sum();              // (n_cells, 1)
            // auto dm_dt_capture_MNPG = MNPG_MNP_ratio.colwise() * dm_dt_capture_MNP_sum;  // (n_cells, n_MNPG_components)

            auto m_from_MNP_sum = m_l_MNP.rowwise().sum() * mask_capture_uncapture +
                                  m_l_MNP.rowwise().sum() * (1 - mask_capture_uncapture);  // (n_cells, 1)
            auto m_from_MNPG = m_l_MNPG.colwise() * mask_capture_uncapture +
                               m_sl_MNPG.colwise() * (1 - mask_capture_uncapture);  // (n_cells, n_MNPG_components)

            auto MNPG_MNP_ratio = m_from_MNPG.colwise() / (m_from_MNP_sum + 1e-20);  // (n_cells, n_MNPG_components)

            auto dm_dt_capture_MNP_sum = dm_dt_capture_MNP.rowwise().sum();              // (n_cells, 1)
            auto dm_dt_capture_MNPG = MNPG_MNP_ratio.colwise() * dm_dt_capture_MNP_sum;  // (n_cells, n_MNPG_components)

            dm_dt_l_MNP -= dm_dt_capture_MNP;
            dm_dt_sl_MNP += dm_dt_capture_MNP;

            dm_dt_l_NMC -= dm_dt_capture_NMC;
            dm_dt_sl_NMC += dm_dt_capture_NMC;

            dm_dt_l_MNPG -= dm_dt_capture_MNPG;
            dm_dt_sl_MNPG += dm_dt_capture_MNPG;

#if LOG_ENABLED
            if (call_count < LOG_FIRST_N_CALLS || call_count % LOG_EVERY_N_CALLS == 0) {
                std::ostringstream oss;

                auto porosity = V_l.array() / cell_volume_l_and_sl;

                oss << "dV_dt_sl = \n" << dV_dt_sl.transpose() << "\n";
                oss << "dm_dt_capture_NMC = \n" << dm_dt_capture_NMC << "\n";
                oss << "dm_dt_capture_MNPG = \n" << dm_dt_capture_MNPG << "\n";
                LOG("rhs_capture.log", oss.str());
            }
#endif
        }  // if (enable_otherPhysics)
    });  // BENCHMARK(t_capture)

    // Convection in liquid phase
    realtype t_convection = 0.0;
    BENCHMARK(t_convection, {
        if (enable_convection) {
            // TODO: WENO still applicable which changing fluid volume between cells -> try out
            if (flowRate > 0.0) {
                internalNonUniformUpwind(m_l, dm_dt_l, flowRate, V_l);
            }
        }
    });

    // Dispersion in liquid phase
    realtype t_dispersion = 0.0;
    BENCHMARK(t_dispersion, {
        if (enable_dispersion) {
            if (dispersion_coefficient > 0.0) {
                // internalUniformDispersion(m_l, dm_dt_l, dispersion_coefficient, cell_distance);
                internalNonUniformDispersion(m_l, dm_dt_l, V_l, dispersion_coefficient, cell_distance);
            }
        }
    });

    realtype t_reactions = 0.0;
    BENCHMARK(t_reactions, {
        if (enable_reactions) {
            // Reactions in liquid phase
            c_l = m_l.colwise() / V_l;

            c_l = c_l.cwiseMax(0);  // Clamp negatives to 0
            // c = c.cwiseAbs();     // Take absolute values

            activities_l.setZero();
            dc_dt_l.setZero();

            reactionSystem.rhs(t, c_l, activities_l, dc_dt_l);
            dm_dt_l += dc_dt_l.colwise() * V_l;

            // Reactions in slurry phase
            c_sl = m_sl.colwise() / V_sl;

            const auto small_slurry_volume_mask = (V_sl.array() / cell_volume_l_and_sl < 0.01)
                                                      .cast<realtype>()
                                                      .rowwise()
                                                      .replicate(n_components());  // (n_cells, 1)

            // If slurry volume is very small, set concentrations to zero to avoid numerical issues and ignore reactions there
            c_sl = small_slurry_volume_mask.select(Array::Zero(n_cells_per_phase, n_components()), c_sl);
            // c_sl = (!c_sl.isFinite()).select(0.0, c_sl);

            c_sl = c_sl.cwiseMax(0);  // Clamp negatives to 0
            // c = c.cwiseAbs();     // Take absolute values

            activities_sl.setZero();
            dc_dt_sl.setZero();

            reactionSystem.rhs(t, c_sl, activities_sl, dc_dt_sl);
            dm_dt_sl += dc_dt_sl.colwise() * V_sl;

            // Diffusion between liquid and slurry NMC
            // internalUniformDispersion(m_l_NMC, dm_dt_l_NMC, dispersion_coefficient, cell_distance);

            // Set fluid concentrations to zero where slurry volume is very small to avoid numerical issues and ignore
            // diffusion there c_l = small_slurry_volume_mask.select(Array::Zero(n_cells_per_phase, n_components()), c_l);

            const auto c_l_NMC = c_l.middleCols(NMC_range.first, NMC_range.second - NMC_range.first + 1);
            const auto c_sl_NMC = c_sl.middleCols(NMC_range.first, NMC_range.second - NMC_range.first + 1);
            const auto contact_area_l_sl = (2 * n_disks * crossSectionArea) / n_cells_per_phase;
            const auto thickness_sl = V_sl / contact_area_l_sl;  // (n_cells, 1)
            if (V_sl.rows() != n_cells_per_phase) {
                throw std::runtime_error("V_sl has incorrect number of rows!");
            }
            auto dm_dt_diffusion_sl_to_l_NMC_raw = dispersion_coefficient_sl_to_l_NMC * contact_area_l_sl *
                                                   ((c_sl_NMC - c_l_NMC).colwise() *
                                                    (thickness_sl).cwiseInverse());  // (n_cells, n_NMC_components)

            const auto dm_dt_diffusion_sl_to_l_NMC = small_slurry_volume_mask
                                                         .middleCols(NMC_range.first, NMC_range.second - NMC_range.first + 1)
                                                         .select(Array::Zero(n_cells_per_phase,
                                                                             NMC_range.second - NMC_range.first + 1),
                                                                 dm_dt_diffusion_sl_to_l_NMC_raw);

            dm_dt_l_NMC += dm_dt_diffusion_sl_to_l_NMC;
            dm_dt_sl_NMC -= dm_dt_diffusion_sl_to_l_NMC;

#if LOG_ENABLED
            if (call_count < LOG_FIRST_N_CALLS || call_count % LOG_EVERY_N_CALLS == 0) {
                std::ostringstream oss;
                oss << "c_l_NMC = \n" << c_l_NMC << "\n";
                oss << "c_sl_NMC = \n" << c_sl_NMC << "\n";
                oss << "dm_dt_diffusion_sl_to_l_NMC = \n" << dm_dt_diffusion_sl_to_l_NMC << "\n\n";
                LOG("rhs_capture.log", oss.str());
            }
#endif
        }
    });

#if BENCHMARK_ENABLED
    {
        std::ostringstream oss;
        oss.setf(std::ios::scientific);
        oss.precision(6);
        oss << "t=" << t << " t_capture=" << t_capture << " t_convection=" << t_convection
            << " t_dispersion=" << t_dispersion << " t_reactions=" << t_reactions << "\n";
        LOG_BENCHMARK("RS_MagneticCaptureProcessChamber_rhs.log", oss.str());
    }
#endif
#if LOG_ENABLED
    // Log non-finite values in c_l
    if (!c_l.allFinite() || (c_l.array() < 0).any()) {
        std::ostringstream oss;
        oss << "Non-finite or negative values in c_l in RS_MagneticProcessChamber at time " << t << ":\n";
        oss << c_l << "\n";
        LOG("unallowed_values.log", oss.str());
    }

    // Log non-finite values in c_sl
    if (!c_sl.allFinite() || (c_sl.array() < 0).any()) {
        std::ostringstream oss;
        oss << "Non-finite or negative values in c_sl in RS_MagneticProcessChamber at time " << t << ":\n";
        oss << c_sl << "\n";
        LOG("unallowed_values.log", oss.str());
    }

    if ((V_l.array() <= 0).any()) {
        std::ostringstream oss;
        oss << "Warning: Negative or zero value detected in V_l in RS_MagneticProcessChamber at time t = " << t
            << ". V_l = " << V_l.transpose() << "\n";
        LOG("unallowed_values.log", oss.str());
    }
    call_count++;
#endif
}