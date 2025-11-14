#ifndef MASS_ACTION_LAW_HPP
#define MASS_ACTION_LAW_HPP

#include "Components/ComponentSystem.hpp"
#include "Reactions/Reaction.hpp"
#include "Reactions/ReactionSystem.hpp"

/** \brief Standard mass action law RHS builder (forward/backward rates). */
std::function<void(realtype t, const ConstArrayMap& concentrations, const ConstArrayMap& activities, ArrayMap& dc_dt)> massActionLaw_rhs(
    std::vector<size_t> reactant_terms,
    std::vector<size_t> product_terms,
    realtype k_forward,
    realtype k_backward);

/** \brief Damped mass action law RHS builder (rate clipping/smoothing via tau). */
std::function<void(realtype t, const ConstArrayMap& concentrations, const ConstArrayMap& activities, ArrayMap& dc_dt)> massActionLaw_damped_rhs(
    std::vector<size_t> reactant_terms,
    std::vector<size_t> product_terms,
    realtype k_forward,
    realtype k_backward,
    realtype tau_reaction);

/** \brief Inverse rate prediction RHS builder (equilibrium extent / tau heuristic). */
std::function<void(realtype t, const ConstArrayMap& concentrations, const ConstArrayMap& activities, ArrayMap& dc_dt)> massActionLaw_inverseRatePrediction_rhs(
    std::vector<size_t> reactant_terms,
    std::vector<size_t> product_terms,
    realtype K_eq,
    realtype tau_reaction);

/** \brief Error function measuring pC deviation (e.g. pH error surrogate). */
std::function<realtype(realtype t, const ConstArrayMap& concentrations, const ConstArrayMap& activities)> massActionLaw_pC_errorFunction(
    size_t idx_H_plus,
    std::vector<size_t> reactant_terms,
    std::vector<size_t> product_terms,
    realtype K_eq);

std::pair<std::vector<size_t>, std::vector<size_t>> parse_massActionLaw(const ComponentSystem& componentSystem,
                                                                        const std::string& reactionStr);
Reaction massActionLaw(const ComponentSystem& componentSystem,
                       const std::string& reactionStr,
                       realtype k_forward,
                       realtype k_backward,
                       int error_component_idx = -1);  // -1 means no error function

Reaction massActionLaw_damped(const ComponentSystem& componentSystem,
                              const std::string& reactionStr,
                              realtype k_forward,
                              realtype k_backward,
                              realtype tau_reaction,
                              int error_component_idx = -1);  // -1 means no error function

Reaction massActionLaw_inverseRatePrediction(const ComponentSystem& componentSystem,
                                             const std::string& reactionStr,
                                             realtype K_eq,
                                             realtype tau_reaction,
                                             int error_component_idx = -1);  // -1 means no error function

#endif  // MASS_ACTION_LAW_HPP