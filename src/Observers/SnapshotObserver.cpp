#include "Observers/SnapshotObserver.hpp"

#include <cstddef>
#include <iostream>

#include "Solver.hpp"
#include "UnitOperations/UnitOperationBase.hpp"
#include "cnpy.h"

void SnapshotObserver::write(realtype t, const realtype* y) {
    t_measured = t;

    // Explicitly copy the data from the mapped view to ensure ownership
    auto mapped_view = mapper(y);
    snapshot = Array(mapped_view);  // Force copy constructor

    LOG("snapshot_observer.log",
        "Observer snapshot at t=" << t_desired << " Snapshot written with size: " << snapshot.size()
                                  << " (rows: " << snapshot.rows() << ", cols: " << snapshot.cols() << ")\n");

    if (compute_errors) {
        error = mapper.unitOperation_ptr->errorFunction(t, y);
    }
}
