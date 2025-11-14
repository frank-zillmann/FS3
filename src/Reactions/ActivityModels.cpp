#include "Reactions/ActivityModels.hpp"

#include <cmath>
#include <string>

#include "Logger.hpp"
#include "PhysicalConstants.hpp"

void NoActivityModel::operator()(realtype t, const ConstArrayMap& concentrations, ArrayMap& activities) const {
    // Vectorized assignment: activities = concentrations for all cells and components at once
    activities = concentrations;
#if LOG_ENABLED
    LOG("activity_models.log",
        "NoActivityModel called at t=" + std::to_string(t) + ": activities set equal to concentrations.\n");
#endif
}

TruesdellJonesActivityModel::TruesdellJonesActivityModel(const ComponentSystem& componentSystem)
    : ActivityModelBase(componentSystem) {
    // B = std::sqrt((2 * constants::elementary_charge * constants::elementary_charge * constants::avogadro) /
    //               (componentSystem.relative_permittivity * constants::vacuum_permittivity * constants::boltzmann *
    //                componentSystem.temperature)) *
    //     std::sqrt(componentSystem.density);

    // A = (constants::elementary_charge * constants::elementary_charge * B) /
    //     (std::log(10) * 8 * constants::pi * constants::vacuum_permittivity * componentSystem.relative_permittivity *
    //      constants::boltzmann * componentSystem.temperature) *
    //     std::sqrt(componentSystem.density);

    // LOG("activity_models.log", "Truesdell-Jones parameters initialized: A = " + std::to_string(A) +
    //                                " [m^1.5/mol^0.5], B = " + std::to_string(B) + " [(m^1.5)/(mol^0.5*m)]\n");

    // A = 1.82483e6 * std::sqrt(componentSystem.density) /
    //     std::pow(componentSystem.relative_permittivity * constants::vacuum_permittivity * componentSystem.temperature, 1.5);
    // A = A * std::pow(10, 1.5);

    // B = 50.2916 * std::sqrt(componentSystem.density / (componentSystem.relative_permittivity *
    //                                                    constants::vacuum_permittivity * componentSystem.temperature));
    // B = B * std::pow(10, 1.5) * 1e9;

    // LOG("activity_models.log", "Truesdell-Jones parameters initialized: A = " + std::to_string(A) +
    //                                " [m^1.5/mol^0.5], B = " + std::to_string(B) + " [(m^1.5)/(mol^0.5*m)]\n");

    // Note: Still some uncertainties / inconsistencies with the units and formulas for A and B in the literature.

    A = 0.5085 * std::pow(10, -1.5);       // Convert from (l^1.5/mol^0.5) to (m^1.5/mol^0.5)
    B = 3.281 * std::pow(10, -1.5) * 1e9;  // Convert from (l/mol)^0.5 to (m^3/mol)^0.5 and from 1/nm to 1/m

    // Attention: Wrong units in paper!
    LOG("activity_models.log", "Truesdell-Jones parameters initialized: A = " + std::to_string(A) +
                                   " [m^1.5/mol^0.5], B = " + std::to_string(B) + " [(m^1.5)/(mol^0.5*m)]\n");
}

void TruesdellJonesActivityModel::operator()(realtype t, const ConstArrayMap& concentrations, ArrayMap& activities) const {
    // Concentrations and Activities may be the same object, so implementation needs to be careful and consider that
    const auto& components = componentSystem.components;
    const auto& n_components = componentSystem.n_components;
    const auto n_cells = concentrations.rows();

    // Calculate ionic strength for each cell (vectorized)
    ColVector ionic_strength = ColVector::Zero(n_cells);
    for (sunindextype componentIndex = 0; componentIndex < n_components; ++componentIndex) {
        if (concentrations.col(componentIndex).hasNaN()) {
            continue;
        }
        const auto& charge = components[componentIndex].charge;
        ionic_strength += 0.5 * concentrations.col(componentIndex) * (charge * charge);
    }

    auto ionic_strength_sqrt = ionic_strength.sqrt();

    // Calculate activity coefficients for each component and cell
    for (sunindextype componentIndex = 0; componentIndex < n_components; ++componentIndex) {
        if (!components[componentIndex].truesdell_jones_alpha.has_value() ||
            !components[componentIndex].truesdell_jones_beta.has_value()) {
            // No activity coefficient correction
            activities.col(componentIndex) = concentrations.col(componentIndex);
            continue;
        }

        const auto& alpha_i = components[componentIndex].truesdell_jones_alpha.value();
        const auto& beta_i = components[componentIndex].truesdell_jones_beta.value();
        const auto& z_i = components[componentIndex].charge;

        // Vectorized calculation of activity coefficients across all cells
        auto numerator = -A * z_i * z_i * ionic_strength_sqrt;
        auto denominator = 1.0 + B * alpha_i * ionic_strength_sqrt;
        auto beta_term = beta_i * ionic_strength;
        auto lg_gamma = numerator / denominator + beta_term;
        auto gamma = (lg_gamma * std::log(10.0)).exp();
#if LOG_ENABLED
        if (call_count < LOG_FIRST_N_CALLS || call_count % LOG_EVERY_N_CALLS == 0) {
            LOG("activity_models.log",
                std::to_string(call_count) + ". Call at t=" + std::to_string(t) + " Component " +
                    components[componentIndex].name + ": ionic_strength mean=" + std::to_string(ionic_strength.mean()) +
                    ", numerator mean=" + std::to_string(numerator.mean()) + ", denominator mean=" +
                    std::to_string(denominator.mean()) + ", beta_term mean=" + std::to_string(beta_term.mean()) +
                    ", lg_gamma mean=" + std::to_string(lg_gamma.mean()) + ", gamma mean=" + std::to_string(gamma.mean()) +
                    " (min: " + std::to_string(gamma.minCoeff()) + ", max: " + std::to_string(gamma.maxCoeff()) + ")\n");
        }
#endif
        activities.col(componentIndex) = gamma * concentrations.col(componentIndex);
    }
#if LOG_ENABLED
    call_count++;
#endif
}