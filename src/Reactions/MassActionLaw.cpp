#include "Reactions/MassActionLaw.hpp"

#include <sundials/sundials_types.h>

#include <cassert>
#include <cmath>
#include <csignal>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <regex>
#include <string>
#include <unordered_map>
#include <utility>

#include "Components/ComponentSystem.hpp"
#include "EigenDataTypes.hpp"
#include "Logger.hpp"
#include "Reactions/Reaction.hpp"

std::function<void(realtype t, const ConstArrayMap& concentrations, const ConstArrayMap& activities, ArrayMap& dc_dt)> massActionLaw_rhs(
    std::vector<size_t> reactant_terms,
    std::vector<size_t> product_terms,
    realtype k_forward,
    realtype k_backward) {
    return [reactant_terms, product_terms, k_forward, k_backward,
            call_count = 0](realtype t, const ConstArrayMap& concentrations, const ConstArrayMap& activities,
                            ArrayMap& dc_dt) mutable {
        const auto n_cells = activities.rows();
        const auto n_components = activities.cols();

        // Calculate forward and backward rates for all cells simultaneously
        ColVector rate_fwd = ColVector::Constant(n_cells, k_forward);
        ColVector rate_bwd = ColVector::Constant(n_cells, k_backward);

        // Calculate forward reaction rates (vectorized)
        // Each index appears as many times as its stoichiometric coefficient
        for (const auto& idx : reactant_terms) {
            rate_fwd *= activities.col(idx);
        }

        // Calculate backward reaction rates (vectorized)
        // Each index appears as many times as its stoichiometric coefficient
        for (const auto& idx : product_terms) {
            rate_bwd *= activities.col(idx);
        }

        auto net_rate = rate_fwd - rate_bwd;

        // Apply stoichiometric changes
        for (size_t idx : reactant_terms) {
            dc_dt.col(idx) -= net_rate;
        }
        for (size_t idx : product_terms) {
            dc_dt.col(idx) += net_rate;
        }

#if LOG_ENABLED
        ++call_count;
        if (call_count < LOG_FIRST_N_CALLS || call_count % LOG_EVERY_N_CALLS == 0) {
            std::ostringstream oss;
            oss.precision(6);
            oss << std::scientific;
            oss << "Mass action law reaction call " << call_count << " at t = " << t << " seconds\n";
            oss << "Reactant indices: ";
            for (size_t idx : reactant_terms) oss << idx << " ";
            oss << "\n";
            oss << "Product indices: ";
            for (size_t idx : product_terms) oss << idx << " ";
            oss << "\n";
            oss << "Forward rate constants: " << k_forward << ", Backward rate constants: " << k_backward << "\n";
            oss << "Activities:\n" << activities << "\n";
            oss << "Forward rate:\n" << rate_fwd.transpose() << "\n";
            oss << "Backward rate:\n" << rate_bwd.transpose() << "\n";
            oss << "Net rate:\n" << net_rate.transpose() << "\n";
            oss << "\n\n";

            LOG("massActionLaw_rhs.log", oss.str());
        }
#endif
    };
}

std::function<void(realtype t, const ConstArrayMap& concentrations, const ConstArrayMap& activities, ArrayMap& dc_dt)> massActionLaw_damped_rhs(
    std::vector<size_t> reactant_terms,
    std::vector<size_t> product_terms,
    realtype k_forward,
    realtype k_backward,
    realtype tau_reaction) {
    return [reactant_terms, product_terms, k_forward, k_backward, tau_reaction,
            call_count = 0](realtype t, const ConstArrayMap& concentrations, const ConstArrayMap& activities,
                            ArrayMap& dc_dt) mutable {
        const auto n_cells = activities.rows();
        const auto n_components = activities.cols();

        // Calculate forward and backward rates for all cells simultaneously
        ColVector rate_fwd = ColVector::Constant(n_cells, k_forward);
        ColVector rate_bwd = ColVector::Constant(n_cells, k_backward);
        ColVector avg_concentration = ColVector::Constant(n_cells, 1);

        ColVector jacobian_sum = ColVector::Constant(n_cells, 0);
        ColVector jacobian_term = ColVector::Constant(n_cells, 0);

        // Calculate forward reaction rates (vectorized)
        // Each index appears as many times as its stoichiometric coefficient
        for (const auto& idx : reactant_terms) {
            rate_fwd *= activities.col(idx);

            avg_concentration = avg_concentration * concentrations.col(idx);

            jacobian_term.setConstant(k_forward);
            int power = 0;
            for (const auto& other_idx : reactant_terms) {
                if (other_idx != idx) {
                    jacobian_term *= activities.col(other_idx);
                }
                if (other_idx == idx) {
                    power++;
                }
            }
            assert(power >= 1 && "Power of component should be at least 1");
            jacobian_term *= power;
            for (int i = 1; i < power; ++i) {
                jacobian_term *= activities.col(idx);
            }
            jacobian_sum += jacobian_term;
        }

        // Calculate backward reaction rates (vectorized)
        // Each index appears as many times as its stoichiometric coefficient
        for (const auto& idx : product_terms) {
            rate_bwd *= activities.col(idx);

            avg_concentration = avg_concentration * concentrations.col(idx);

            jacobian_term.setConstant(k_backward);
            int power = 0;
            for (const auto& other_idx : product_terms) {
                if (other_idx != idx) {
                    jacobian_term *= activities.col(other_idx);
                }
                if (other_idx == idx) {
                    power++;
                }
            }
            assert(power >= 1 && "Power of component should be at least 1");
            jacobian_term *= power;
            for (int i = 1; i < power; ++i) {
                jacobian_term *= activities.col(idx);
            }
            jacobian_sum += jacobian_term;
        }

        std::size_t n_reactant_and_products = reactant_terms.size() + product_terms.size();
        realtype root = 1.0 / static_cast<realtype>(n_reactant_and_products);
        avg_concentration = avg_concentration.pow(root);

        auto net_rate = rate_fwd - rate_bwd;

        // auto c_magnitude = 0.2e-3;
        // realtype max_flowRate = 1.79583e-05;
        // realtype min_cell_volume = 7.38299e-06;
        // realtype alpha_reaction = 1e-3;
        // auto max_reaction_rate = alpha_reaction * (max_flowRate / min_cell_volume) * c_magnitude;

        auto max_rate = avg_concentration / tau_reaction;

        // auto max_rate = ColVector::Constant(n_cells, tau_reaction);

        // Damping to avoid numerical instabilities due to very high reaction rates

        // Bound the net rate to be within [-max_rate, max_rate]
        // auto bounded_net_rate = net_rate.cwiseMax(-max_rate).cwiseMin(max_rate);
        // auto bounded_net_rate = (net_rate >= 0).select((net_rate + 1).log(), -((-net_rate + 1).log()));
        // auto bounded_net_rate = (tau_reaction / jacobian_sum) * net_rate;

        // auto max_rate = tau_reaction;

        auto bounded_net_rate = net_rate / (1.0 + (net_rate.cwiseAbs() / max_rate));
        // auto bounded_net_rate = (net_rate.array() > 0)
        //                             .select(ColVector::Constant(n_cells, max_rate),
        //                                     (net_rate.array() < 0)
        //                                         .select(ColVector::Constant(n_cells, -max_rate),
        //                                                 ColVector::Constant(n_cells, 0)));

        // Apply stoichiometric changes
        for (size_t idx : reactant_terms) {
            dc_dt.col(idx) -= bounded_net_rate;
        }
        for (size_t idx : product_terms) {
            dc_dt.col(idx) += bounded_net_rate;
        }

        // TODO: Try damping with excluding H2O

#if LOG_ENABLED
        ++call_count;
        if (call_count < LOG_FIRST_N_CALLS || call_count % LOG_EVERY_N_CALLS == 0) {
            std::ostringstream oss;
            oss.precision(6);
            oss << std::scientific;
            oss << "Mass action law reaction call " << call_count << " at t = " << t << " seconds\n";
            oss << "Reactant indices: ";
            for (size_t idx : reactant_terms) oss << idx << " ";
            oss << "\n";
            oss << "Product indices: ";
            for (size_t idx : product_terms) oss << idx << " ";
            oss << "\n";
            oss << "Forward rate constants: " << k_forward << ", Backward rate constants: " << k_backward
                << ", Damping factor: " << tau_reaction << "\n";
            oss << "Activities:\n" << activities << "\n";
            oss << "Forward rate:\n" << rate_fwd.transpose() << "\n";
            oss << "Backward rate:\n" << rate_bwd.transpose() << "\n";
            oss << "Net rate:\n" << net_rate.transpose() << "\n";
            // oss << "Jacobian sum:\n" << jacobian_sum.transpose() << "\n";
            oss << "Average concentration:\n" << avg_concentration.transpose() << "\n";
            oss << "Bounded net rate:\n" << bounded_net_rate.transpose() << "\n";
            oss << "\n\n";

            LOG("massActionLaw_damped_rhs.log", oss.str());
        }
#endif
    };
}

std::function<void(realtype t, const ConstArrayMap& concentrations, const ConstArrayMap& activities, ArrayMap& dc_dt)> massActionLaw_inverseRatePrediction_rhs(
    std::vector<size_t> reactant_terms,
    std::vector<size_t> product_terms,
    realtype K_eq,
    realtype tau_reaction) {
    return [reactant_terms, product_terms, K_eq, tau_reaction,
            call_count = 0](realtype t, const ConstArrayMap& concentrations, const ConstArrayMap& activities,
                            ArrayMap& dc_dt) mutable {
        const auto n_cells = activities.rows();
        const auto n_components = activities.cols();

        if (reactant_terms.size() == 1 && product_terms.size() == 2) {
            // A -> B + C
            size_t idx_A = reactant_terms[0];
            size_t idx_B = product_terms[0];
            size_t idx_C = product_terms[1];

            // k_f * (A - r) = k_b * (B + r)(C + r)
            // k_f * A - k_f * r = k_b * r^2 + k_b * (B + C) * r + k_b * B * C
            // r^2 + (B + C + K_eq) * r + (B * C - K_eq * A) = 0 with K_eq = k_f / k_b

            auto p = activities.col(idx_B) + activities.col(idx_C) + K_eq;
            auto q = activities.col(idx_B) * activities.col(idx_C) - K_eq * activities.col(idx_A);

            auto radicant = (p.square() - 4.0 * q).sqrt();

            const auto delta_xi_1 = (-p + radicant) / 2.0;

            // Assuming delta_xi_1 is the physically correct solution
            auto net_rate = delta_xi_1 / tau_reaction;
            dc_dt.col(idx_A) -= net_rate;
            dc_dt.col(idx_B) += net_rate;
            dc_dt.col(idx_C) += net_rate;
#if LOG_ENABLED || !defined(NDEBUG)
            // const auto residual_current = K_eq * activities.col(idx_A) - activities.col(idx_B) * activities.col(idx_C);
            const auto Q_current = activities.col(idx_B) * activities.col(idx_C) / activities.col(idx_A);
            const auto Q_relative_error_current_raw = (Q_current - K_eq).cwiseAbs() / K_eq;

            const auto c_A_eq_1 = activities.col(idx_A) - delta_xi_1;
            const auto c_B_eq_1 = activities.col(idx_B) + delta_xi_1;
            const auto c_C_eq_1 = activities.col(idx_C) + delta_xi_1;
            // const auto residual_eq_1 = k_forward * c_A_eq_1 - k_backward * c_B_eq_1 * c_C_eq_1;
            const auto Q_eq_1 = c_B_eq_1 * c_C_eq_1 / c_A_eq_1;
            // Set NaN/Inf values to 0 to avoid assertion failures when dividing by zero
            const auto Q_relative_error_1_raw = (Q_eq_1 - K_eq).cwiseAbs() / K_eq;
            const auto Q_relative_error_1 = Q_relative_error_1_raw.isFinite().select(Q_relative_error_1_raw, 0.0);

            const auto delta_xi_2 = (-p - radicant) / 2.0;
            const auto c_A_eq_2 = activities.col(idx_A) - delta_xi_2;
            const auto c_B_eq_2 = activities.col(idx_B) + delta_xi_2;
            const auto c_C_eq_2 = activities.col(idx_C) + delta_xi_2;
            // const auto residual_eq_2 = k_forward * c_A_eq_2 - k_backward * c_B_eq_2 * c_C_eq_2;
            const auto Q_eq_2 = c_B_eq_2 * c_C_eq_2 / c_A_eq_2;
            // Set NaN/Inf values to 0 to avoid assertion failures when dividing by zero
            const auto Q_relative_error_2_raw = (Q_eq_2 - K_eq).cwiseAbs() / K_eq;
            const auto Q_relative_error_2 = Q_relative_error_2_raw.isFinite().select(Q_relative_error_2_raw, 0.0);
#endif
#if LOG_ENABLED
            if (call_count < LOG_FIRST_N_CALLS || call_count % LOG_EVERY_N_CALLS == 0) {
                std::ostringstream oss;
                oss.precision(6);
                oss << std::scientific;
                oss << "Mass action law reaction call " << call_count << " at t = " << t << " seconds\n";
                oss << "Reactant index: " << idx_A << ", Product indices: " << idx_B << ", " << idx_C << "\n";
                oss << "K_eq: " << K_eq << "\n";
                oss << "Activities:\n" << activities << "\n";
                // oss << "residual_current:\n" << residual_current.transpose() << "\n";
                oss << "Q_current:\n" << Q_current.transpose() << "\n";
                oss << "Q_relative_error_current:\n" << Q_relative_error_current_raw.transpose() << "\n";
                oss << "p:\n" << p.transpose() << "\n";
                oss << "q:\n" << q.transpose() << "\n";
                oss << "Radicant:\n" << radicant.transpose() << "\n";
                oss << "delta_xi_1:\n" << delta_xi_1.transpose() << "\n";
                oss << "A_eq_1:\n" << c_A_eq_1.transpose() << "\n";
                oss << "B_eq_1:\n" << c_B_eq_1.transpose() << "\n";
                oss << "C_eq_1:\n" << c_C_eq_1.transpose() << "\n";
                // oss << "residual_eq_1:\n" << residual_eq_1.transpose() << "\n";
                oss << "Q_eq_1:\n" << Q_eq_1.transpose() << "\n";
                oss << "Relative error Q_eq_1 vs K_eq:\n" << Q_relative_error_1.transpose() << "\n";
                oss << "delta_xi_2:\n" << delta_xi_2.transpose() << "\n";
                oss << "A_eq_2:\n" << c_A_eq_2.transpose() << "\n";
                oss << "B_eq_2:\n" << c_B_eq_2.transpose() << "\n";
                oss << "C_eq_2:\n" << c_C_eq_2.transpose() << "\n";
                // oss << "residual_eq_2:\n" << residual_eq_2.transpose() << "\n";
                oss << "Q_eq_2:\n" << Q_eq_2.transpose() << "\n";
                oss << "Relative error Q_eq_2 vs K_eq:\n" << Q_relative_error_2.transpose() << "\n";
                oss << "dr/dt:\n" << net_rate.transpose() << "\n\n";

                LOG("massActionLaw_inverseRatePrediction_rhs.log", oss.str());
            }
#endif
            assert(((c_A_eq_1 + 1e-12) >= 0).all() && ((c_B_eq_1 + 1e-12) >= 0).all() &&
                   ((c_C_eq_1 + 1e-12) >= 0).all() && "Equilibrium concentrations (solution 1) must be non-negative");

            assert(((c_A_eq_2 - 1e-12) < 0).all() || ((c_B_eq_2 - 1e-12) < 0).all() ||
                   ((c_C_eq_2 - 1e-12) < 0).all() &&
                       "Equilibrium concentrations (solution 2) should contain negative values");

            // Skip these assertions because they are fragile e.g. for extremely dilute solutions, etc.

            // assert((Q_relative_error_1 < 1e-1).all() &&
            //        "Equilibrium concentrations (solution 1) must satisfy equilibrium condition");

            // assert((residual_eq_1.abs() < 1e-6).all() &&
            //        "Equilibrium concentrations (solution 1) must satisfy equilibrium condition");  // Quite high numerical errors due to extreme K_eq values

            // assert((residual_eq_1.abs() < residual_current.abs()).all() &&
            //        "Equilibrium concentrations (solution 1) must be closer to equilibrium than current concentrations");

        } else {
            throw std::runtime_error("Mass action law reaction currently only supports 1 reactant to 2 products.");
        }
        ++call_count;
    };
}

std::function<realtype(realtype t, const ConstArrayMap& concentrations, const ConstArrayMap& activities)> massActionLaw_pC_errorFunction(
    size_t idx_H_plus,
    std::vector<size_t> reactant_terms,
    std::vector<size_t> product_terms,
    realtype K_eq) {
    return [idx_H_plus, reactant_terms, product_terms, K_eq,
            call_count = 0](realtype t, const ConstArrayMap& concentrations,
                            const ConstArrayMap& activities) mutable -> realtype {
        // TODO: If H+ is a reactant, swap reactants and products and invert K_eq

        const auto n_cells = activities.rows();

        ColVector max_concentration = ColVector::Constant(n_cells, 0);
        ColVector c_H_plus_acitivity_approx = ColVector::Constant(n_cells, K_eq);

        // Each index appears as many times as its stoichiometric coefficient
        for (size_t idx : reactant_terms) {
            c_H_plus_acitivity_approx *= activities.col(idx);
            max_concentration = max_concentration.cwiseMax(concentrations.col(idx));
        }

        // Calculate backward reaction rates (vectorized)
        // Each index appears as many times as its stoichiometric coefficient
        for (size_t idx : product_terms) {
            if (idx != idx_H_plus) {
                c_H_plus_acitivity_approx /= activities.col(idx);
            }
            max_concentration = max_concentration.cwiseMax(concentrations.col(idx));
        }

        auto pH_approx = -c_H_plus_acitivity_approx.log10() + 3;  // Adjust for concentrations being mol/m³ instead of mol/L
        auto pH = -activities.col(idx_H_plus).log10() + 3;  // Adjust for concentrations being mol/m³ instead of mol/L
        auto pH_error = (pH - pH_approx).eval();

        // Only consider pH_error for cells where max_concentration > 1e-12, to avoid high errors in pure water / extremely dilute solutions
        auto mask = (max_concentration.array() > 1e-12);
        pH_error = mask.select(pH_error, 0);
        // Set NaN values to 0 because they arise if other products and reactants are simultaneously 0, so the reaction is irrelevant
        pH_error = pH_error.unaryExpr([](realtype v) { return std::isfinite(v) ? v : 0.0; });

#if LOG_ENABLED
        if (call_count < LOG_FIRST_N_CALLS || call_count % LOG_EVERY_N_CALLS == 0) {
            std::ostringstream oss;
            oss.precision(6);
            oss << std::scientific;

            oss << "Mass action law reaction call " << call_count << " at t = " << t << " seconds\n";
            oss << "Reactant indices: ";
            for (size_t idx : reactant_terms) oss << idx << " ";
            oss << "\n";
            oss << "Product indices: ";
            for (size_t idx : product_terms) oss << idx << " ";
            oss << "\n";

            oss << "Equilibrium constant K_eq: " << K_eq << "\n";
            oss << "Activities:\n" << activities << "\n";

            oss << "pH:" << pH.transpose() << "\n";
            oss << "pH_approx: " << pH_approx.transpose() << "\n";
            oss << "pH_error: " << pH_error.transpose() << "\n\n";
            LOG("massActionLaw_pC_errorFunction.log", oss.str());
        }
        call_count++;
#endif

        return pH_error.maxCoeff();
    };
}

// Helper to trim whitespace
static std::string trim(const std::string& s) {
    auto start = s.find_first_not_of(" \t");
    auto end = s.find_last_not_of(" \t");
    return (start == std::string::npos) ? "" : s.substr(start, end - start + 1);
}

// Parse a side of the reaction (e.g., "2*A + B") into a map of component name -> stoichiometry
static std::unordered_map<std::string, int> parseSide(const std::string& side) {
    std::unordered_map<std::string, int> result;

    std::regex term_regex(R"((\d*)\s*\*?\s*([A-Za-z_][A-Za-z0-9_⁺⁻₀₁₂₃₄₅₆₇₈₉+-]*))");
    std::regex delimiter_regex(R"(\s+\+\s+)");  // Matches ' + '

    std::sregex_token_iterator iter(side.begin(), side.end(), delimiter_regex, -1);
    std::sregex_token_iterator end;
    std::vector<std::string> terms(iter, end);

    for (auto& term : terms) {
        term = trim(term);
        std::smatch match;
        if (std::regex_match(term, match, term_regex)) {
            int coeff = match[1].str().empty() ? 1 : std::stoi(match[1].str());
            std::string name = match[2];
            result[name] += coeff;
        }
    }
    return result;
}

std::pair<std::vector<size_t>, std::vector<size_t>> parse_massActionLaw(
    const ComponentSystem& componentSystem,
    const std::string& reactionStr) {  // Split reaction string into left and right sides
    std::string delimiter = "<=>";
    auto pos = reactionStr.find(delimiter);
    if (pos == std::string::npos) {
        throw std::invalid_argument("Reaction string must contain '<=>' delimiter.");
    }
    std::string lhs = trim(reactionStr.substr(0, pos));
    std::string rhs = trim(reactionStr.substr(pos + delimiter.size()));

    auto reactants = parseSide(lhs);
    auto products = parseSide(rhs);

    // Prepare vectors of indices where each index appears as many times as its stoichiometric coefficient
    std::vector<size_t> reactant_terms, product_terms;
    for (const auto& [name, coeff] : reactants) {
        size_t idx = componentSystem.getIdx(name);
        for (int i = 0; i < coeff; ++i) {
            reactant_terms.push_back(idx);
        }
    }
    for (const auto& [name, coeff] : products) {
        size_t idx = componentSystem.getIdx(name);
        for (int i = 0; i < coeff; ++i) {
            product_terms.push_back(idx);
        }
    }

#if LOG_ENABLED
    std::ostringstream oss;
    oss << "Reaction string: " << reactionStr << "\n";
    oss << "Parsed reactants:\n";
    for (const auto& [name, coeff] : reactants) {
        oss << "  " << name << ": " << coeff << "\n";
    }
    oss << "Parsed products:\n";
    for (const auto& [name, coeff] : products) {
        oss << "  " << name << ": " << coeff << "\n";
    }
    oss << "Reactant indices: ";
    for (size_t idx : reactant_terms) oss << idx << " ";
    oss << "\n";
    oss << "Product indices: ";
    for (size_t idx : product_terms) oss << idx << " ";
    oss << "\n\n";
    LOG("parse_massActionLaw.log", oss.str());
#endif

    return std::make_pair(reactant_terms, product_terms);
}

Reaction massActionLaw(const ComponentSystem& componentSystem,
                       const std::string& reactionStr,
                       realtype k_forward,
                       realtype k_backward,
                       int error_component_idx) {
    auto [reactant_terms, product_terms] = parse_massActionLaw(componentSystem, reactionStr);

    if (error_component_idx >= 0) {
        // Calculate equilibrium constant K_eq
        auto K_eq = k_forward / k_backward;

        return Reaction(massActionLaw_rhs(reactant_terms, product_terms, k_forward, k_backward),
                        massActionLaw_pC_errorFunction(error_component_idx, reactant_terms, product_terms, K_eq));
    } else {
        return Reaction(massActionLaw_rhs(reactant_terms, product_terms, k_forward, k_backward));
    }
}

Reaction massActionLaw_damped(const ComponentSystem& componentSystem,
                              const std::string& reactionStr,
                              realtype k_forward,
                              realtype k_backward,
                              realtype tau_reaction,
                              int error_component_idx) {
    auto [reactant_terms, product_terms] = parse_massActionLaw(componentSystem, reactionStr);

    if (error_component_idx >= 0) {
        // Calculate equilibrium constant K_eq
        auto K_eq = k_forward / k_backward;

        return Reaction(massActionLaw_damped_rhs(reactant_terms, product_terms, k_forward, k_backward, tau_reaction),
                        massActionLaw_pC_errorFunction(error_component_idx, reactant_terms, product_terms, K_eq));
    } else {
        return Reaction(massActionLaw_damped_rhs(reactant_terms, product_terms, k_forward, k_backward, tau_reaction));
    }
}

Reaction massActionLaw_inverseRatePrediction(const ComponentSystem& componentSystem,
                                             const std::string& reactionStr,
                                             realtype K_eq,
                                             realtype tau_reaction,
                                             int error_component_idx) {
    auto [reactant_terms, product_terms] = parse_massActionLaw(componentSystem, reactionStr);

    if (error_component_idx >= 0) {
        // Calculate equilibrium constant K_eq

        return Reaction(massActionLaw_inverseRatePrediction_rhs(reactant_terms, product_terms, K_eq, tau_reaction),
                        massActionLaw_pC_errorFunction(error_component_idx, reactant_terms, product_terms, K_eq));
    } else {
        return Reaction(massActionLaw_inverseRatePrediction_rhs(reactant_terms, product_terms, K_eq, tau_reaction));
    }
}
