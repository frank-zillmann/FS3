#include "Observers/PHConverter.hpp"

#include <cmath>
#include <stdexcept>

#include "Components/ComponentSystem.hpp"
#include "EigenDataTypes.hpp"
#include "Observers/TimeSeriesObserver.hpp"

ColVector convert_to_pH(const TimeSeriesObserver& observer,
                        std::size_t cell_index,
                        realtype cell_volume,
                        const ComponentSystem& componentSystem,
                        const ActivityModelBase& activityModel) {
    // Get all snapshots as a vector (need to cast away const as method is not const)
    auto [obs_vector, shape] = const_cast<TimeSeriesObserver&>(observer).get_all_snapshots_as_vector();
    auto [n_timesteps, n_cells, n_components] = shape;

    // Validate inputs
    if (obs_vector.empty()) {
        throw std::runtime_error("TimeSeriesObserver contains no data");
    }
    if (cell_index >= n_cells) {
        throw std::runtime_error("cell_index is out of bounds for the observer data");
    }
    if (cell_volume <= 0) {
        throw std::runtime_error("cell_volume must be positive");
    }

    // Validate component system matches observer data
    if (componentSystem.n_components != n_components) {
        throw std::runtime_error("Number of components in observer does not match component system");
    }

    // Find the index of H⁺ component
    std::size_t idx_H_plus = componentSystem.getIdx("H⁺");

    // Allocate arrays for processing
    // Shape: (n_timesteps, n_components)
    Array masses(n_timesteps, n_components);
    Array concentrations(n_timesteps, n_components);
    Array activities(n_timesteps, n_components);

    // Extract data for the specified cell from all timesteps
    for (std::size_t t = 0; t < n_timesteps; ++t) {
        // Calculate the starting index for this timestep and cell
        std::size_t base_idx = t * n_cells * n_components + cell_index * n_components;

        // Copy the component data for this cell at this timestep
        for (std::size_t comp = 0; comp < n_components; ++comp) {
            masses(t, comp) = obs_vector[base_idx + comp];
        }
    }

    // Convert masses (kg / mol) to concentrations (kg/m³ / mol/m³)
    concentrations = masses / cell_volume;

    // Create maps for activity model (it expects mapped arrays)
    // Note: We use stride (n_components, 1) since data is row-major contiguous
    PhaseStride stride(n_components, 1);
    ConstArrayMap concentrations_map(concentrations.data(), n_timesteps, n_components, stride);
    ArrayMap activities_map(activities.data(), n_timesteps, n_components, stride);

    // Compute activities using the activity model
    activityModel(0.0, concentrations_map, activities_map);

    // Calculate pH values from H⁺ activity
    // pH = -log10(a_H⁺) where a_H⁺ is in mol/L
    // Since our concentrations are in mol/m³, we need to add 3 to convert:
    // pH = -log10(a_H⁺ [mol/m³]) + 3 = -log10(a_H⁺ [mol/L])
    ColVector pH_values = -(activities.col(idx_H_plus).log10()) + 3.0;

    return pH_values;
}
