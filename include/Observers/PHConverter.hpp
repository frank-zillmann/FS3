#ifndef PH_CONVERTER_HPP
#define PH_CONVERTER_HPP

#include <cstddef>

#include "Components/ComponentSystem.hpp"
#include "EigenDataTypes.hpp"
#include "Observers/TimeSeriesObserver.hpp"
#include "Reactions/ActivityModels.hpp"
#include "sundials/sundials_types.h"

/**
 * @brief Convert TimeSeriesObserver data to pH values for a specific cell
 *
 * This function extracts mass data from a TimeSeriesObserver for a specific cell,
 * converts masses to molar concentrations, computes activities using the provided
 * activity model, and calculates pH values based on H⁺ activity.
 *
 * @param observer TimeSeriesObserver containing mass data over time
 * @param cell_index Index of the cell to extract data from
 * @param cell_volume Volume of the cell in m³
 * @param componentSystem Component system to convert masses to concentrations
 * @param activityModel Activity model to compute activities from concentrations
 * @return ColVector Array of pH values with shape (n_timesteps,)
 */
ColVector convert_to_pH(const TimeSeriesObserver& observer,
                        std::size_t cell_index,
                        realtype cell_volume,
                        const ComponentSystem& componentSystem,
                        const ActivityModelBase& activityModel);

#endif  // PH_CONVERTER_HPP
