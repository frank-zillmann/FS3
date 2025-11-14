#include <sundials/sundials_types.h>

#include "Components/ComponentSystem.hpp"
#include "EigenDataTypes.hpp"
#include "Logger.hpp"
#include "Reactions/ActivityModels.hpp"
#include "Reactions/Reaction.hpp"

template <EigenArrayLike CArrayType, WritableEigenArrayLike ActivityArrayType, WritableEigenArrayLike DcDtArrayType>
void ReactionSystem::rhs(realtype t, const CArrayType& c, ActivityArrayType& activities, DcDtArrayType& dc_dt) const {
    const auto n_cells = c.rows();
    const auto n_components = c.cols();

    // Set negative concentrations to zero and log warnings
#if LOG_ENABLED
    for (Eigen::Index cell = 0; cell < n_cells; ++cell) {
        for (Eigen::Index comp = 0; comp < n_components; ++comp) {
            if (c(cell, comp) < 0) {
                const auto& component = componentSystem(comp);
                LOG("reaction_system_negativeConcentration.log",
                    "Negative concentration for component " + component.name + " in cell " + std::to_string(cell) +
                        " at time " + std::to_string(t) + ": " + std::to_string(c(cell, comp)) + "\n");
            }
        }
    }
#endif

    // Maps for access for activity model and reactions
    ConstArrayMap c_map(c.data(), n_cells, n_components, PhaseStride(n_components, 1));
    ArrayMap activities_map(activities.data(), n_cells, n_components,
                            PhaseStride(n_components, 1));  // For activity models that modify activities
    ConstArrayMap const_activities_map(activities.data(), n_cells, n_components,
                                       PhaseStride(n_components, 1));  // For reactions that do not modify activities
    ArrayMap dc_dt_map(dc_dt.data(), n_cells, n_components, PhaseStride(n_components, 1));

    // Call activity model to calculate activities based on concentrations
    activityModel(t, c_map, activities_map);

    // Apply all reactions
    int reaction_count = 0;
    for (const auto& reaction : reactions) {
        reaction.rhs(t, c_map, const_activities_map, dc_dt_map);

#if LOG_ENABLED
        realtype error = reaction.errorFunction(t, c_map, const_activities_map);
        LOG("reactionSystem_rhs.log", "t = " + std::to_string(t) + "\t\tReaction No. = " + std::to_string(reaction_count) +
                                          ":\t\terror = " + std::to_string(error) + +"\n");
#endif
        reaction_count++;
    }
}

template <EigenArrayLike CArrayType, WritableEigenArrayLike ActivityArrayType>
realtype ReactionSystem::errorFunction(realtype t, const CArrayType& c, ActivityArrayType& activities) const {
    const auto n_cells = c.rows();
    const auto n_components = c.cols();

    // Set negative concentrations to zero and log warnings
#if LOG_ENABLED
    for (Eigen::Index cell = 0; cell < n_cells; ++cell) {
        for (Eigen::Index comp = 0; comp < n_components; ++comp) {
            if (c(cell, comp) < 0) {
                const auto& component = componentSystem(comp);
                LOG("reaction_system_negativeConcentration.log",
                    "Negative concentration for component " + component.name + " in cell " + std::to_string(cell) +
                        " at time " + std::to_string(t) + ": " + std::to_string(c(cell, comp)) + "\n");
            }
        }
    }
#endif

    // Maps for access for activity model and reactions
    ConstArrayMap c_map(c.data(), n_cells, n_components, PhaseStride(n_components, 1));
    ArrayMap activities_map(activities.data(), n_cells, n_components,
                            PhaseStride(n_components, 1));  // For activity models that modify activities
    ConstArrayMap const_activities_map(activities.data(), n_cells, n_components,
                                       PhaseStride(n_components, 1));  // For reactions that do not modify activities

    // Call activity model to calculate activities based on concentrations
    activityModel(t, c_map, activities_map);

    realtype current_max_error = 0.0;
    // Apply all reactions
    int reaction_count = 0;
    for (const auto& reaction : reactions) {
        realtype error = reaction.errorFunction(t, c_map, const_activities_map);

        if (error > max_error) {
            LOG("reactionSystem_errorFunction.log",
                "Reaction error exceeded previous max at time " + std::to_string(t) + ": reaction " +
                    std::to_string(reaction_count) + " error = " + std::to_string(error) + "\n");
            max_error = error;
        }

        LOG("reactionSystem_errorFunction.log", "t = " + std::to_string(t) +
                                                    "\t\tReaction No. = " + std::to_string(reaction_count) +
                                                    ":\t\tcurrent_error = " + std::to_string(error) + +"\n");

        if (error > current_max_error) {
            current_max_error = error;
        }

        ++reaction_count;
    }

    return current_max_error;
}