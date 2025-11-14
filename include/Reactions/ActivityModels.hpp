#ifndef ACTIVITY_MODELS_HPP
#define ACTIVITY_MODELS_HPP

#include <sundials/sundials_types.h>

#include <cassert>
#include <cmath>

#include "Components/ComponentSystem.hpp"
#include "EigenDataTypes.hpp"

/**
 * @brief Base class for activity coefficient models
 *
 * Computes activities from concentrations, vectorized across cells.
 */
class ActivityModelBase {
   public:
    ActivityModelBase(const ComponentSystem& componentSystem) : componentSystem(componentSystem) {}
    virtual ~ActivityModelBase() = default;

    // Eigen-based interface for vectorized operations across multiple cells
    virtual void operator()(realtype t, const ConstArrayMap& concentrations, ArrayMap& activities) const = 0;

   protected:
    const ComponentSystem& componentSystem;
};

/**
 * @brief Identity activity model (activities = concentrations)
 */
class NoActivityModel : public ActivityModelBase {
   public:
    NoActivityModel(const ComponentSystem& componentSystem) : ActivityModelBase(componentSystem) {}
    virtual void operator()(realtype t, const ConstArrayMap& concentrations, ArrayMap& activities) const override;
};

/**
 * @brief Truesdell–Jones activity model for electrolytes
 *
 * Extended Debye–Hückel form with finite ion size and specific interaction terms.
 */
class TruesdellJonesActivityModel : public ActivityModelBase {
   public:
    TruesdellJonesActivityModel(const ComponentSystem& componentSystem);

    virtual void operator()(realtype t, const ConstArrayMap& concentrations, ArrayMap& activities) const override;

   protected:
    realtype A;
    realtype B;

#if LOG_ENABLED
    mutable int call_count = 0;
#endif
};

#endif  // ACTIVITY_MODELS_HPP