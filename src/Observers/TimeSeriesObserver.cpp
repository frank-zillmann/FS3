#include "Observers/TimeSeriesObserver.hpp"

#include <cstddef>
#include <iostream>

#include "Solver.hpp"
#include "UnitOperations/UnitOperationBase.hpp"
#include "cnpy.h"

TimeSeriesObserver::TimeSeriesObserver(const realtype t_start,
                                       const realtype t_end,
                                       const std::size_t n_snapshots,
                                       const ArrayMapper& mapper,
                                       bool compute_errors) {
    if (n_snapshots < 2) {
        throw std::runtime_error("TimeSeriesObserver must contain 2 or more snapshots.");
    }
    snapshotObserver_ptrs.reserve(n_snapshots);

    realtype step = (t_end - t_start) / static_cast<realtype>(n_snapshots - 1);
    for (std::size_t i = 0; i < n_snapshots - 1; ++i) {
        realtype t = t_start + step * static_cast<realtype>(i);
        snapshotObserver_ptrs.emplace_back(std::make_shared<SnapshotObserver>(t, mapper, compute_errors));
    }
    snapshotObserver_ptrs.emplace_back(std::make_shared<SnapshotObserver>(t_end, mapper, compute_errors));
}

void TimeSeriesObserver::add_all_snapshots_to_solver(Solver& solver) {
    for (const auto& snapshotObserver_ptr : snapshotObserver_ptrs) {
        solver.add(snapshotObserver_ptr);
    }
}

std::pair<std::vector<realtype>, std::array<size_t, 3>> TimeSeriesObserver::get_all_snapshots_as_vector() {
    const std::size_t n_snapshotObservers = snapshotObserver_ptrs.size();
    const std::size_t n_cells = snapshotObserver_ptrs[0]->snapshot.rows();
    const std::size_t n_components = snapshotObserver_ptrs[0]->snapshot.cols();

    std::size_t n_time_steps = 0;
    for (const auto& snapshotObserver_ptr : snapshotObserver_ptrs) {
        if (snapshotObserver_ptr->snapshot.rows() == n_cells && snapshotObserver_ptr->snapshot.cols() == n_components) {
            n_time_steps++;
        } else if (snapshotObserver_ptr->snapshot.size() == 0) {
            // Skip empty snapshots
        } else {
            throw std::runtime_error("Inconsistent snapshot sizes in TimeSeriesObserver.");
        }
    }
    if (n_time_steps == 0) {
        std::cerr << "Warning: No valid snapshots in TimeSeriesObserver." << std::endl;
        return {{}, {0, 0, 0}};
    }
    if (n_time_steps != n_snapshotObservers) {
        std::cerr << "Warning: Some snapshots were skipped due to inconsistent sizes in TimeSeriesObserver." << std::endl;
    }

    std::vector<realtype> obs_vector(n_time_steps * n_cells * n_components);

    for (size_t t_idx = 0; t_idx < n_time_steps; ++t_idx) {
        const auto& snapshot = snapshotObserver_ptrs[t_idx]->snapshot;
        std::memcpy(obs_vector.data() + (t_idx * n_cells * n_components), snapshot.data(),
                    sizeof(realtype) * n_cells * n_components);
    }

    return {obs_vector, {n_time_steps, n_cells, n_components}};
}

void TimeSeriesObserver::save_to_npz(std::string zipname, std::string fname, std::string mode) {
    auto [obs_vector, shape] = get_all_snapshots_as_vector();
    auto [n_time_steps, n_cells, n_components] = shape;

    if (obs_vector.empty()) {
        std::cerr << "Warning: No data to save in TimeSeriesObserver." << std::endl;
        return;
    }

    cnpy::npz_save(zipname, fname, obs_vector.data(),
                   {static_cast<size_t>(n_time_steps), static_cast<size_t>(n_cells), static_cast<size_t>(n_components)},
                   mode);

    LOG("time_series_observer.log", "Saved TimeSeriesObserver to " << zipname << " with variable name " << fname
                                                                   << " having shape (" << n_time_steps << ", "
                                                                   << n_cells << ", " << n_components << ")\n");
}

void TimeSeriesObserver::save_to_npy(std::string filename) {
    auto [obs_vector, shape] = get_all_snapshots_as_vector();
    auto [n_time_steps, n_cells, n_components] = shape;

    if (obs_vector.empty()) {
        std::cerr << "Warning: No data to save in TimeSeriesObserver." << std::endl;
        return;
    }

    cnpy::npy_save(filename, obs_vector.data(),
                   {static_cast<size_t>(n_time_steps), static_cast<size_t>(n_cells), static_cast<size_t>(n_components)});

    LOG("time_series_observer.log", "Saved TimeSeriesObserver to " << filename << " having shape (" << n_time_steps
                                                                   << ", " << n_cells << ", " << n_components << ")\n");
}