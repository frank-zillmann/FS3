#ifndef TIMESERIESOBSERVER_HPP
#define TIMESERIESOBSERVER_HPP

#include <cstddef>
#include <cstring>
#include <memory>
#include <vector>

#include "EigenDataTypes.hpp"
#include "Observers/SnapshotObserver.hpp"
#include "sundials/sundials_types.h"

class TimeSeriesObserver {
   public:
    TimeSeriesObserver(const realtype t_start,
                       const realtype t_end,
                       const std::size_t n_snapshots,
                       const ArrayMapper& mapper,
                       bool compute_errors = false);

    TimeSeriesObserver(const realtype t_start,
                       const realtype t_end,
                       const std::size_t n_snapshots,
                       const ArrayMapper& mapper,
                       Solver& solver,
                       bool add_to_solver = true,
                       bool compute_errors = false)
        : TimeSeriesObserver(t_start, t_end, n_snapshots, mapper, compute_errors) {
        if (add_to_solver) {
            add_all_snapshots_to_solver(solver);
        }
    }

    void add_all_snapshots_to_solver(Solver& solver);

    // Returns: [obs_vector, [n_time_steps, n_cells, n_components]]
    std::pair<std::vector<realtype>, std::array<size_t, 3>> get_all_snapshots_as_vector();

    void save_to_npz(std::string zipname, std::string fname, std::string mode);

    void save_to_npy(std::string filename);

   private:
    std::vector<std::shared_ptr<SnapshotObserver>> snapshotObserver_ptrs;
};

#endif  // TIMESERIESOBSERVER_HPP