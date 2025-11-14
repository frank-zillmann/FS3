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
    realtype t_desired;
    realtype t_measured;
    ArrayMapper mapper;
    Array snapshot;
    bool compute_errors;
    realtype error;

    SnapshotObserver(const realtype time, const ArrayMapper& mapper, bool compute_errors = false)
        : t_desired(time), mapper(mapper), compute_errors(compute_errors), error(-1.0) {
        snapshot.resize(mapper.rows, mapper.cols);
    }

    void write(realtype t, const realtype* y);
};

#endif  // SNAPSHOTOBSERVER_HPP