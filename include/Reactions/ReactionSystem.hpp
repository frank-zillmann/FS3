#ifndef REACTION_SYSTEM_HPP
#define REACTION_SYSTEM_HPP

#include <vector>

#include "EigenDataTypes.hpp"
#include "sundials/sundials_types.h"

#include "Reactions/Reaction.hpp"

// Forward declarations
class ComponentSystem;
class ActivityModelBase;

/**
 * @brief Collection of reactions with an activity model
 *
 * Applies the activity model, evaluates all reaction RHS functors, and
 * accumulates dc/dt and optional error metrics for diagnostics.
 */
class ReactionSystem {
   public:
    ReactionSystem(const ComponentSystem& componentSystem, const ActivityModelBase& activityModel)
        : componentSystem(componentSystem), activityModel(activityModel) {}

    void add(const Reaction& reaction);

    // Expects concentrations in SI units (mol/m^3 or kg/m^3)
    template <EigenArrayLike CArrayType, WritableEigenArrayLike ActivityArrayType, WritableEigenArrayLike DcDtArrayType>
    void rhs(realtype t, const CArrayType& c, ActivityArrayType& activities, DcDtArrayType& dc_dt) const;

    // Expects concentrations in SI units (mol/m^3 or kg/m^3)
    template <EigenArrayLike CArrayType, WritableEigenArrayLike ActivityArrayType>
    realtype errorFunction(realtype t, const CArrayType& c, ActivityArrayType& activities) const;

    std::vector<Reaction> reactions;
    const ComponentSystem& componentSystem;
    const ActivityModelBase& activityModel;

    mutable realtype max_error = 0.0;  // Maximum error encountered during simulation
};

#include "Reactions/ReactionSystem.tpp"

#endif  // REACTION_SYSTEM_HPP