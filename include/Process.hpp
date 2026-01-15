#ifndef PROCESS_HPP
#define PROCESS_HPP

#include <sundials/sundials_types.h>

#include <memory>
#include <vector>

#include "Components/ComponentSystem.hpp"
#include "EigenDataTypes.hpp"
#include "UnitOperations/UnitOperationBase.hpp"
class Solver;

/**
 * @brief Aggregates unit operations and flow connections
 *
 * Holds ordered unit operations, connection descriptors (source/target mappers
 * with time-dependent flow rate) and provides process-level RHS assembly helpers.
 */
class Process {
    friend class Solver;

   public:
    Process(const ComponentSystem& componentSystem,
            const std::initializer_list<std::shared_ptr<UnitOperationBase>> unitOperations,
            realtype t_end = std::numeric_limits<realtype>::infinity())
        : componentSystem(componentSystem), unitOperations(unitOperations), t_end(t_end) {
        c_inlet_buffer.resize(1, componentSystem.n_components);
        dc_dt_inlet_buffer.resize(1, componentSystem.n_components);
        dc_dt_outlet_buffer.resize(1, componentSystem.n_components);
    }

    void addConnection(const ArrayMapper from, const ArrayMapper to, std::function<realtype(realtype)> flowRateFunction) {
        connections.push_back({from, to, flowRateFunction});
    }

    void rhs_connections(realtype t, const realtype* y, realtype* dy_dt, const Solver& solver) const;

    const std::vector<std::shared_ptr<UnitOperationBase>>& getUnitOperations() const { return unitOperations; }

   protected:
    // List of all components in the process in the same order as in y
    // -> maps relative index to detailed information about the component
    const ComponentSystem& componentSystem;

    // List of all unit operations in the process
    const std::vector<std::shared_ptr<UnitOperationBase>> unitOperations;

    realtype t_end;  // End time of the process (default: infinity)

    struct Connection {
        ArrayMapper fromMapper;  // helper to get the data from the from unit operation
        ArrayMapper toMapper;    // helper to get the data to the to unit operation
        std::function<realtype(realtype)> flowRateFunction;
    };
    std::vector<Connection> connections;  // List of all connections between unit operations in the process

    mutable RowVector c_inlet_buffer;
    mutable RowVector dc_dt_inlet_buffer;
    mutable RowVector dc_dt_outlet_buffer;
};

#endif  // PROCESS_HPP