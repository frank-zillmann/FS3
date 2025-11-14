#ifndef CONVECTION_DISPERSION_OPERATORS_TPP
#define CONVECTION_DISPERSION_OPERATORS_TPP

#include <cassert>
#include <sstream>
#include <stdexcept>

#include "EigenDataTypes.hpp"
#include "Logger.hpp"

template <EigenArrayLike MArrayType, WritableEigenArrayLike DmDtArrayType>
void internalUniformUpwind(const MArrayType& m, DmDtArrayType& dm_dt, realtype flowRate, realtype cell_volume) {
    const auto nCells = m.rows();
    if (nCells < 2) {
        return;  // No upwind needed for single cell
    }

    if (flowRate < 0.0) {
        throw std::runtime_error("Flow rate must not be negative for upwind convection.");
    }

    auto flux = (flowRate / cell_volume) * m.topRows(nCells - 1);  // m(0…nCells-2, :) / cell_volume for concentration

    dm_dt.topRows(nCells - 1) -= flux;
    dm_dt.bottomRows(nCells - 1) += flux;

#if LOG_ENABLED
    static size_t call_count = 0;
    if (call_count < LOG_FIRST_N_CALLS || call_count % LOG_EVERY_N_CALLS == 0) {
        std::ostringstream oss;
        oss << call_count << ": internalUniformUpwind\n";
        oss << "flowRate = " << flowRate << ", volume(cell_volume) = " << cell_volume << "\n";
        oss << "m(topRows) = \n" << m.topRows(nCells - 1) << "\n";
        oss << "flux = \n" << flux << "\n";
        LOG("internalUniformUpwind.log", oss.str());
    }
    ++call_count;
#endif
}

template <EigenArrayLike MArrayType, WritableEigenArrayLike DmDtArrayType, EigenArrayLike V_l_ArrayType>
void internalNonUniformUpwind(const MArrayType& m, DmDtArrayType& dm_dt, realtype flowRate, const V_l_ArrayType& V_l) {
    const auto nCells = m.rows();
    if (nCells < 2) {
        return;  // No upwind needed for single cell
    }

    if (flowRate < 0.0) {
        throw std::runtime_error("Flow rate must not be negative for upwind convection.");
    }

    if (V_l.cols() != 1 || V_l.rows() != nCells) {
        throw std::runtime_error("Cell volumes array must be a column vector with size equal to number of cells.");
    }

    auto flux = m.topRows(nCells - 1).colwise() * (flowRate * V_l.head(nCells - 1).cwiseInverse());
    // Alternative: m.topRows(nCells - 1).array() / V_l.head(nCells - 1);

    dm_dt.topRows(nCells - 1) -= flux;
    dm_dt.bottomRows(nCells - 1) += flux;

#if LOG_ENABLED
    static size_t call_count = 0;
    if (call_count < LOG_FIRST_N_CALLS || call_count % LOG_EVERY_N_CALLS == 0) {
        std::ostringstream oss;
        oss << call_count << ": internalNonUniformUpwind\n";
        oss << "flowRate = " << flowRate << "\n";
        oss << "V_l(head) = " << V_l.head(nCells - 1).transpose() << "\n";
        oss << "m(topRows) = \n" << m.topRows(nCells - 1) << "\n";
        oss << "flux = \n" << flux << "\n";
        LOG("internalNonUniformUpwind.log", oss.str());
    }
    ++call_count;
#endif
}

template <EigenArrayLike MArrayType, WritableEigenArrayLike DmDtArrayType>
void internalUniformWENO3(const MArrayType& m, DmDtArrayType& dm_dt, realtype flowRate, realtype cell_volume) {
    const auto nCells = m.rows();
    if (nCells < 3) {
        throw std::runtime_error("WENO3 requires at least 3 cells for reconstruction.");
    }
    if (flowRate < 0.0) {
        throw std::runtime_error("Flow rate must not be negative for WENO3 convection.");
    }

    // WENO3 coefficients
    constexpr realtype eps = 1e-6;
    // Stencil: c_{i-1}, c_i, c_{i+1}

    // Convert masses to concentrations for WENO reconstruction
    auto c_im1 = m.topRows(nCells - 2) / cell_volume;
    auto c_i = m.middleRows(1, nCells - 2) / cell_volume;
    auto c_ip1 = m.bottomRows(nCells - 2) / cell_volume;

    // Ensure dimensions match in Debug mode
    assert(c_im1.cols() == c_i.cols() && c_i.cols() == c_ip1.cols());
    assert(c_im1.rows() == nCells - 2 && c_i.rows() == nCells - 2 && c_ip1.rows() == nCells - 2);

    // Compute candidate fluxes (linear weights)
    auto q0 = 0.5 * (c_ip1 + c_i);
    auto q1 = 0.5 * (-c_im1 + 3 * c_i);

    // Smoothness indicators
    auto beta0 = (c_ip1 - c_i).square();
    auto beta1 = (c_i - c_im1).square();

    // Linear weights
    constexpr realtype gamma0 = 1.0 / 3.0;
    constexpr realtype gamma1 = 2.0 / 3.0;

    // Nonlinear weights
    auto alpha0 = gamma0 / ((eps + beta0).square());
    auto alpha1 = gamma1 / ((eps + beta1).square());
    auto alpha_sum = alpha0 + alpha1;
    auto w0 = alpha0 / alpha_sum;
    auto w1 = alpha1 / alpha_sum;

    // WENO flux at i+1/2 (concentration flux, then multiply by flowRate for mass flux)
    auto flux = flowRate * (w0 * q0 + w1 * q1);

    dm_dt.middleRows(1, nCells - 2) -= flux;
    dm_dt.bottomRows(nCells - 2) += flux;

    // Upwind between cell 0 and cell 1 (inlined calculation)
    auto flux_0to1 = (flowRate / cell_volume) * m.row(0);
    dm_dt.row(0) -= flux_0to1;
    dm_dt.row(1) += flux_0to1;

#if LOG_ENABLED
    static size_t call_count = 0;
    if (call_count < LOG_FIRST_N_CALLS || call_count % LOG_EVERY_N_CALLS == 0) {
        std::ostringstream oss;
        oss << call_count << ": internalUniformWENO3\n";
        oss << "flowRate = " << flowRate << ", volume(cell_volume) = " << cell_volume << "\n";
        oss << "m = \n" << m << "\n";
        oss << "flux(mid) = \n" << flux << "\n";
        oss << "flux(0->1) = \n" << flux_0to1 << "\n";
        LOG("internalUniformWENO3.log", oss.str());
    }
    ++call_count;
#endif
}

template <EigenArrayLike MFromArrayType, WritableEigenArrayLike DmDtFromArrayType, EigenArrayLike MToArrayType, WritableEigenArrayLike DmDtToArrayType>
void connectionUpwind(const MFromArrayType& m_from,
                      DmDtFromArrayType& dm_dt_from,
                      realtype cell_volume_from,
                      const MToArrayType& m_to,
                      DmDtToArrayType& dm_dt_to,
                      realtype flowRate) {
    if (m_from.rows() != 1 || m_to.rows() != 1 || dm_dt_from.rows() != 1 || dm_dt_to.rows() != 1) {
        throw std::runtime_error("Connection upwind requires 1D arrays for m_from, dm_dt_from, m_to, and dm_dt_to.");
    }

    auto flux = (flowRate / cell_volume_from) * m_from.row(0);

    dm_dt_from.row(0) -= flux;
    dm_dt_to.row(0) += flux;

#if LOG_ENABLED
    static size_t call_count = 0;
    if (call_count < LOG_FIRST_N_CALLS || call_count % LOG_EVERY_N_CALLS == 0) {
        std::ostringstream oss;
        oss << call_count << ": connectionUpwind\n";
        oss << "flowRate = " << flowRate << ", volume(cell_volume_from) = " << cell_volume_from << "\n";
        oss << "m_from = " << m_from << "\n";
        oss << "flux = " << flux << "\n";
        LOG("connectionUpwind.log", oss.str());
    }
    ++call_count;
#endif
}

template <EigenArrayLike MArrayType, WritableEigenArrayLike DmDtArrayType>
void internalUniformDispersion(const MArrayType& m,
                               DmDtArrayType& dm_dt,
                               realtype dispersion_coefficient,
                               realtype cell_distance) {
    const auto nCells = m.rows();
    if (nCells < 2) {
        return;  // No dispersion needed for single cell
    }

    // Calculate the dispersion flux
    auto flux = (dispersion_coefficient / (cell_distance * cell_distance)) *
                (m.topRows(nCells - 1) - m.bottomRows(nCells - 1));

    // Update the rate of change
    dm_dt.topRows(nCells - 1) -= flux;
    dm_dt.bottomRows(nCells - 1) += flux;

#if LOG_ENABLED
    static size_t call_count = 0;
    if (call_count < LOG_FIRST_N_CALLS || call_count % LOG_EVERY_N_CALLS == 0) {
        std::ostringstream oss;
        oss << call_count << ": internalUniformDispersion\n";
        oss << "D = " << dispersion_coefficient << ", dx = " << cell_distance << "\n";
        oss << "m = \n" << m << "\n";
        oss << "flux = \n" << flux << "\n";
        LOG("internalUniformDispersion.log", oss.str());
    }
    ++call_count;
#endif
}

template <EigenArrayLike MArrayType, WritableEigenArrayLike DmDtArrayType, EigenArrayLike V_l_ArrayType>
void internalNonUniformDispersion(const MArrayType& m,
                                  DmDtArrayType& dm_dt,
                                  const V_l_ArrayType& V_l,
                                  realtype dispersion_coefficient,
                                  realtype cell_distance) {
    const auto nCells = m.rows();
    if (nCells < 2) return;  // No dispersion needed for single cell

    const auto inv_dx = 1 / cell_distance;

    const auto c = m.colwise() / V_l;

    // interface "areas" (via volumes): A_interface = (V_i + V_d(i)}) / (2 * Δx)
    const auto A_interface = (V_l.head(nCells - 1) + V_l.tail(nCells - 1)) * (0.5 * inv_dx);

    // flux_{i<->d(i)} = -D * A_interface * (c_d(i) - c_i) / Δx
    const auto factor = (-dispersion_coefficient * inv_dx * A_interface);
    const auto flux = (c.bottomRows(nCells - 1) - c.topRows(nCells - 1)).colwise() * factor;

    // mass balance: dm/dt = flux_{u(i)<->i} - flux_{i<->d(i)}
    dm_dt.topRows(nCells - 1) -= flux;
    dm_dt.bottomRows(nCells - 1) += flux;

#if LOG_ENABLED
    static size_t call_count = 0;
    if (call_count < LOG_FIRST_N_CALLS || call_count % LOG_EVERY_N_CALLS == 0) {
        std::ostringstream oss;
        oss << call_count << ": internalNonUniformDispersion\n";
        oss << "D = " << dispersion_coefficient << ", dx = " << cell_distance << "\n";
        oss << "V_l = " << V_l.transpose() << "\n";
        oss << "m = \n" << m << "\n";
        oss << "flux = \n" << flux << "\n";
        LOG("internalNonUniformDispersion.log", oss.str());
    }
    ++call_count;
#endif
}

#endif  // CONVECTION_DISPERSION_OPERATORS_TPP
