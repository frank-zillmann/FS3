#ifndef SNAPSHOTOBSERVER_HPP
#define SNAPSHOTOBSERVER_HPP

#include <cstddef>
#include <cstring>
#include <memory>
#include <vector>

#include "EigenDataTypes.hpp"
#include "sundials/sundials_types.h"

class Solver;  // Forward declaration
class SnapshotObserver {
   public:
    const realtype t_desired;
    const ArrayMapper mapper;
    const bool compute_errors;

    realtype t_measured;
    Array snapshot;
    realtype error;

    SnapshotObserver(const realtype time, const ArrayMapper& mapper, bool compute_errors = false)
        : t_desired(time), mapper(mapper), compute_errors(compute_errors), error(-1.0) {
        snapshot.resize(mapper.rows, mapper.cols);
    }

    void write(realtype t, const realtype* y);
};

#endif  // SNAPSHOTOBSERVER_HPP