#include "Process.hpp"

#include <sundials/sundials_types.h>

#include <iostream>

#include "ConvectionDispersionOperators.hpp"
#include "EigenDataTypes.hpp"
#include "Solver.hpp"
#include "UnitOperations/Inlet.hpp"
#include "UnitOperations/UnitOperationBase.hpp"

void Process::rhs_connections(realtype t, const realtype* y, realtype* dy_dt, const Solver& solver) const {
    for (const auto& connection : connections) {
        realtype flowRate = connection.flowRateFunction(t);
        if (flowRate == 0.0) {
            continue;  // No flow -> skip
        }

        // Using lambda functions because ArrayMap and ConstArrayMap cannot be copied but need to be constructed outside
        // the if else scopes
        // FROM
        auto c_from = [&] {
            if (connection.fromMapper.unitOperation_ptr->getType() == UnitOperationType::Inlet) {
                auto inlet = static_cast<const Inlet*>(connection.fromMapper.unitOperation_ptr);
                c_inlet_buffer = inlet->getSolution(t);
                return connection.fromMapper(static_cast<const realtype*>(c_inlet_buffer.data()));
            } else {
                assert(connection.fromMapper.unitOperation_ptr->getType() != UnitOperationType::Outlet &&
                       "Outlet cannot be a source in a connection.");
                return connection.fromMapper(y + solver.getUnitOperationStartIdx(connection.fromMapper.unitOperation_ptr));
            }
        }();

        auto dc_dt_from = [&] {
            if (connection.fromMapper.unitOperation_ptr->getType() == UnitOperationType::Inlet) {
                return connection.fromMapper(dc_dt_inlet_buffer.data());
            } else {
                return connection.fromMapper(dy_dt +
                                             solver.getUnitOperationStartIdx(connection.fromMapper.unitOperation_ptr));
            }
        }();

        // TO
        auto c_to = [&] {
            if (connection.toMapper.unitOperation_ptr->getType() == UnitOperationType::Outlet) {
                return connection.toMapper(
                    c_from.data());  // Using from concentrations to avoid dispersion / define boundary condition
            } else {
                assert(connection.toMapper.unitOperation_ptr->getType() != UnitOperationType::Inlet &&
                       "Inlet cannot be a sink in a connection.");
                return connection.toMapper(y + solver.getUnitOperationStartIdx(connection.toMapper.unitOperation_ptr));
            }
        }();

        auto dc_dt_to = [&] {
            if (connection.toMapper.unitOperation_ptr->getType() == UnitOperationType::Outlet) {
                dc_dt_outlet_buffer.resize(1.0, componentSystem.n_components);
                dc_dt_outlet_buffer.setConstant(0.0);
                return connection.toMapper(dc_dt_outlet_buffer.data());
            } else {
                return connection.toMapper(dy_dt + solver.getUnitOperationStartIdx(connection.toMapper.unitOperation_ptr));
            }
        }();

        const auto& V_l_from = connection.fromMapper.unitOperation_ptr->get_V_l();

        connectionUpwind(c_from, dc_dt_from, V_l_from(V_l_from.size() - 1), c_to, dc_dt_to, flowRate);
    }
}

// Old:
// for (const auto& connection : connections) {
//     if (t >= connection.startTime && t <= connection.endTime) {
//         const auto fromType = connection.fromMapper.unitOperation_ptr->getType();
//         const auto toType = connection.toMapper.unitOperation_ptr->getType();

//         if (fromType == UnitOperationType::Unknown || toType == UnitOperationType::Unknown) {
//             throw std::runtime_error("Unknown unit operation type in connection.");
//         } else if (fromType == UnitOperationType::Inlet && toType == UnitOperationType::Outlet) {
//             throw std::runtime_error("Inlet to Outlet connections make no sense.");
//         } else if (fromType == UnitOperationType::Outlet) {
//             throw std::runtime_error("Oulet cannot be a source in a connection.");
//         } else if (toType == UnitOperationType::Inlet) {
//             throw std::runtime_error("Inlet cannot be a sink in a connection.");
//         } else if (fromType == UnitOperationType::Inlet) {
//             auto inlet = static_cast<const Inlet*>(connection.fromMapper.unitOperation_ptr);
//             auto y_from = inlet->getSolution(t);
//             if (y_from.size() != componentSystem.n_components) {
//                 std::cerr << "Warning: Inlet solution vector size (" << y_from.size()
//                           << ") is different than number of components (" << componentSystem.n_components
//                           << "). Padding with zeros." << std::endl;
//                 y_from.resize(componentSystem.n_components, 0.0);
//             }
//             const auto to_idx = solver.getUnitOperationStartIdx(connection.toMapper.unitOperation_ptr) +
//             connection.toIndex; const auto& v_to = connection.toMapper.unitOperation_ptr->get_cell_volume();

//             for (std::size_t comp_idx = 0; comp_idx < componentSystem.n_components; ++comp_idx) {
//                 const auto flux = connection.flowRate * y_from[comp_idx];
//                 dy_dt[to_idx + comp_idx] += flux / v_to;
//             }
//         } else if (toType == UnitOperationType::Outlet) {
//             const auto from_idx = solver.getUnitOperationStartIdx(connection.fromMapper.unitOperation_ptr) +
//             connection.fromIndex; const auto& v_from = connection.fromMapper.unitOperation_ptr->get_cell_volume();

//             for (std::size_t comp_idx = 0; comp_idx < componentSystem.n_components; ++comp_idx) {
//                 const auto flux = connection.flowRate * y[from_idx + comp_idx];
//                 dy_dt[from_idx + comp_idx] -= flux / v_from;
//             }
//         } else {  // For all other unit operations
//             const auto from_idx = solver.getUnitOperationStartIdx(connection.fromMapper.unitOperation_ptr) + connection.fromIndex;
//             const auto to_idx = solver.getUnitOperationStartIdx(connection.toMapper.unitOperation_ptr) + connection.toIndex;

//             const auto& v_from = connection.fromMapper.unitOperation_ptr->get_cell_volume();
//             const auto& v_to = connection.toMapper.unitOperation_ptr->get_cell_volume();

//             for (std::size_t comp_idx = 0; comp_idx < componentSystem.n_components; ++comp_idx) {
//                 const auto flux = connection.flowRate * y[from_idx + comp_idx];
//                 dy_dt[from_idx + comp_idx] -= flux / v_from;
//                 dy_dt[to_idx + comp_idx] += flux / v_to;
//             }
//         }
//     }
// }