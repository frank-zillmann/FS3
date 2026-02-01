/**
 * @file bindings.cpp
 * @brief Python bindings for FS3 library using nanobind
 *
 * This file provides Python bindings for the core FS3 classes:
 * - Components: Component, ComponentSystem
 * - UnitOperations: Inlet, Outlet, Volume, Pipe, RS_MagneticCaptureProcessChamber
 * - Reactions: Reaction, ReactionSystem, MassActionLaw functions, ActivityModels
 * - Process and Solver
 * - Observers: SnapshotObserver, TimeSeriesObserver
 */

#include <nanobind/nanobind.h>
#include <nanobind/eigen/dense.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/function.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "Components/Component.hpp"
#include "Components/ComponentSystem.hpp"
#include "EigenDataTypes.hpp"
#include "Observers/PHConverter.hpp"
#include "Observers/SnapshotObserver.hpp"
#include "Observers/TimeSeriesObserver.hpp"
#include "Process.hpp"
#include "Reactions/ActivityModels.hpp"
#include "Reactions/MassActionLaw.hpp"
#include "Reactions/Reaction.hpp"
#include "Reactions/ReactionSystem.hpp"
#include "Solver.hpp"
#include "UnitOperations/Inlet.hpp"
#include "UnitOperations/Outlet.hpp"
#include "UnitOperations/Pipe.hpp"
#include "UnitOperations/RS_MagneticCaptureProcessChamber.hpp"
#include "UnitOperations/Volume.hpp"

namespace nb = nanobind;
using namespace nb::literals;

NB_MODULE(fs3, m) {
    m.doc() = "FS³ - Fast and Flexible Framework for Simple Simulations of Separation-Processes";

    // ==================== Enums ====================

    nb::enum_<Type>(m, "ComponentType")
        .value("NonMagneticComponent", Type::NonMageneticComponent)
        .value("MagneticNanoParticle", Type::MagneticNanoParticle)
        .value("MagneticNanoParticleGroup", Type::MagneticNanoParticleGroup);

    nb::enum_<Phase>(m, "Phase")
        .value("Liquid", Phase::Liquid)
        .value("SolidOrSlurry", Phase::SolidOrSlurry);

    nb::enum_<UnitOperationType>(m, "UnitOperationType")
        .value("Unknown", UnitOperationType::Unknown)
        .value("Inlet", UnitOperationType::Inlet)
        .value("Outlet", UnitOperationType::Outlet)
        .value("Volume", UnitOperationType::Volume)
        .value("Pipe", UnitOperationType::Pipe)
        .value("RS_MagneticCaptureProcessChamber", UnitOperationType::RS_MagneticCaptureProcessChamber);

    nb::enum_<SolverType>(m, "SolverType")
        .value("BDF", SolverType::BDF)
        .value("ADAMS", SolverType::ADAMS)
        .value("ERK", SolverType::ERK)
        .value("ARK", SolverType::ARK);

    // ==================== Component ====================

    nb::class_<Component>(m, "Component")
        .def(nb::init<const std::string&>(), "name"_a)
        // Keyword-argument friendly constructor: allows setting common properties
        .def("__init__", [](Component& self, const std::string& name, nb::kwargs kwargs) {
            new (&self) Component(name);
            if (kwargs.contains("molar_mass_kg_per_mol"))
                self.setMolarMass(nb::cast<double>(kwargs["molar_mass_kg_per_mol"]));
            if (kwargs.contains("charge"))
                self.setCharge(nb::cast<int>(kwargs["charge"]));
            if (kwargs.contains("type"))
                self.setType(nb::cast<Type>(kwargs["type"]));
            if (kwargs.contains("radius_m"))
                self.setRadius(nb::cast<double>(kwargs["radius_m"]));
            if (kwargs.contains("density_kg_per_m3"))
                self.setDensity(nb::cast<double>(kwargs["density_kg_per_m3"]));
            if (kwargs.contains("magnetic_saturation_A_per_m"))
                self.setMagneticSaturation(nb::cast<double>(kwargs["magnetic_saturation_A_per_m"]));
            if (kwargs.contains("truesdell_jones_alpha_meter") && kwargs.contains("truesdell_jones_beta_cubicmeter_per_mol"))
                self.setTruesdellJonesParameters(nb::cast<double>(kwargs["truesdell_jones_alpha_meter"]),
                                                 nb::cast<double>(kwargs["truesdell_jones_beta_cubicmeter_per_mol"]));
            else if (kwargs.contains("truesdell_jones_alpha_meter") || kwargs.contains("truesdell_jones_beta_cubicmeter_per_mol"))
                throw std::invalid_argument(
                    "Both truesdell_jones_alpha_meter and truesdell_jones_beta_cubicmeter_per_mol must be provided together!");
        }, "name"_a, "kwargs"_a)
        .def_rw("name", &Component::name)
        .def_rw("charge", &Component::charge)
        .def_rw("type", &Component::type)
        .def_rw("molar_mass", &Component::molarMass)
        .def_rw("radius", &Component::radius)
        .def_rw("density", &Component::density)
        .def_rw("magnetic_saturation", &Component::magnetic_saturation)
        .def_rw("truesdell_jones_alpha", &Component::truesdell_jones_alpha)
        .def_rw("truesdell_jones_beta", &Component::truesdell_jones_beta)
        .def("__repr__", [](const Component& c) {
            return "<Component '" + c.name + "' charge=" + std::to_string(c.charge) + ">";
        });

    // ==================== ComponentSystem ====================

    nb::class_<ComponentSystem>(m, "ComponentSystem")
        .def(nb::init<const std::vector<Component>&>(), "components"_a)
        .def_ro("n_components", &ComponentSystem::n_components)
        .def_ro("components", &ComponentSystem::components)
        .def_ro("medium_name", &ComponentSystem::medium_name)
        .def_ro("temperature", &ComponentSystem::temperature)
        .def_ro("density", &ComponentSystem::density)
        .def_ro("dynamic_viscosity", &ComponentSystem::dynamic_viscosity)
        .def_ro("relative_permittivity", &ComponentSystem::relative_permittivity)
        .def("get_idx", &ComponentSystem::getIdx, "name"_a)
        .def("get_molar_masses", &ComponentSystem::getMolarMasses)
        .def("get_inv_molar_masses", &ComponentSystem::getInvMolarMasses)
        .def("get_charges", &ComponentSystem::getCharges)
        .def("__getitem__", nb::overload_cast<const std::string&>(&ComponentSystem::operator(), nb::const_), "name"_a)
        .def("__getitem__", nb::overload_cast<sunindextype>(&ComponentSystem::operator(), nb::const_), "index"_a)
        .def("__len__", [](const ComponentSystem& cs) { return cs.n_components; })
        .def("__repr__", [](const ComponentSystem& cs) {
            return "<ComponentSystem with " + std::to_string(cs.n_components) + " components>";
        });

    // ==================== Activity Models ====================

    nb::class_<ActivityModelBase>(m, "ActivityModelBase");

    nb::class_<NoActivityModel, ActivityModelBase>(m, "NoActivityModel")
        .def(nb::init<const ComponentSystem&>(), "component_system"_a);

    nb::class_<TruesdellJonesActivityModel, ActivityModelBase>(m, "TruesdellJonesActivityModel")
        .def(nb::init<const ComponentSystem&>(), "component_system"_a);

    // ==================== Reaction ====================

    // Note: Reaction uses std::function with ArrayMap types which are complex to bind directly.
    // We expose it but creation should go through factory functions like massActionLaw.
    nb::class_<Reaction>(m, "Reaction")
        .def("__repr__", [](const Reaction&) {
            return "<Reaction>";
        });

    // ==================== ReactionSystem ====================

    nb::class_<ReactionSystem>(m, "ReactionSystem")
        .def(nb::init<const ComponentSystem&, const ActivityModelBase&>(),
             "component_system"_a, "activity_model"_a,
             nb::keep_alive<1, 2>(), nb::keep_alive<1, 3>())
        .def("add", &ReactionSystem::add, "reaction"_a)
        .def_ro("reactions", &ReactionSystem::reactions)
        .def("__len__", [](const ReactionSystem& rs) { return rs.reactions.size(); })
        .def("__repr__", [](const ReactionSystem& rs) {
            return "<ReactionSystem with " + std::to_string(rs.reactions.size()) + " reactions>";
        });

    // ==================== Mass Action Law Factory Functions ====================

    m.def("parse_mass_action_law", &parse_massActionLaw,
          "component_system"_a, "reaction_str"_a,
          "Parse a reaction string like 'A + B <-> C + D' and return (reactant_indices, product_indices)");

    m.def("mass_action_law", &massActionLaw,
          "component_system"_a, "reaction_str"_a, "k_forward"_a, "k_backward"_a, "error_component_idx"_a = -1,
          "Create a mass action law reaction from a string like 'H2O <-> H+ + OH-'");

    m.def("mass_action_law_damped", &massActionLaw_damped,
          "component_system"_a, "reaction_str"_a, "k_forward"_a, "k_backward"_a, "tau_reaction"_a, "error_component_idx"_a = -1,
          "Create a damped mass action law reaction with rate smoothing");

    m.def("mass_action_law_inverse_rate_prediction", &massActionLaw_inverseRatePrediction,
          "component_system"_a, "reaction_str"_a, "K_eq"_a, "tau_reaction"_a, "error_component_idx"_a = -1,
          "Create an inverse rate prediction mass action law reaction");

    // ==================== ArrayMapper ====================

    nb::class_<ArrayMapper>(m, "ArrayMapper")
        .def(nb::init<>())
        .def_ro("rows", &ArrayMapper::rows)
        .def_ro("cols", &ArrayMapper::cols)
        .def_ro("start_idx", &ArrayMapper::startIdx)
        .def_prop_ro("size", [](const ArrayMapper& m) { return m.rows * m.cols; })
        .def("__repr__", [](const ArrayMapper& m) {
            return "<ArrayMapper rows=" + std::to_string(m.rows) +
                   " cols=" + std::to_string(m.cols) +
                   " start_idx=" + std::to_string(m.startIdx) + ">";
        });

    // ==================== UnitOperationBase ====================

    nb::class_<UnitOperationBase>(m, "UnitOperationBase")
        .def("get_type", &UnitOperationBase::getType)
        .def("n_components", &UnitOperationBase::n_components)
        .def("y_size", &UnitOperationBase::y_size)
        .def_ro("n_cells", &UnitOperationBase::n_cells)
        .def_rw("y", &UnitOperationBase::y);

    // ==================== Inlet ====================

    nb::class_<Inlet, UnitOperationBase>(m, "Inlet")
        .def(nb::init<const ReactionSystem&, std::function<RowVector(const realtype&)>>(),
             "reaction_system"_a, "solution_function"_a,
             nb::keep_alive<1, 2>())
        .def("out", &Inlet::out)
        .def("get_solution", &Inlet::getSolution, "time"_a)
        .def("__repr__", [](const Inlet&) {
            return "<Inlet>";
        });

    // ==================== Outlet ====================

    nb::class_<Outlet, UnitOperationBase>(m, "Outlet")
        .def(nb::init<const ReactionSystem&>(), "reaction_system"_a, nb::keep_alive<1, 2>())
        .def("in_", &Outlet::in)  // 'in' is a Python keyword
        .def("__repr__", [](const Outlet&) {
            return "<Outlet>";
        });

    // ==================== Volume ====================

    nb::class_<Volume, UnitOperationBase>(m, "Volume")
        .def(nb::init<const ReactionSystem&, realtype>(),
             "reaction_system"_a, "volume"_a = std::numeric_limits<realtype>::quiet_NaN(),
             nb::keep_alive<1, 2>())
        .def("set_const_initial_concentration", &Volume::setConstInitialConcentration, "concentrations"_a)
        .def("in_", &Volume::in)  // 'in' is a Python keyword
        .def("out", &Volume::out)
        .def("all", &Volume::all)
        .def_ro("cell_volume", &Volume::cell_volume)
        .def("__repr__", [](const Volume& v) {
            return "<Volume cell_volume=" + std::to_string(v.cell_volume) + ">";
        });

    // ==================== Pipe ====================

    nb::class_<Pipe, UnitOperationBase>(m, "Pipe")
        .def(nb::init<const ReactionSystem&, sunindextype, realtype, realtype,
                      std::function<realtype(realtype)>, realtype>(),
             "reaction_system"_a, "n_cells"_a, "cross_section_area"_a, "length"_a,
             "flow_rate_function"_a, "dispersion_coefficient"_a,
             nb::keep_alive<1, 2>())
        .def("set_const_initial_concentration", &Pipe::setConstInitialConcentration, "concentrations"_a)
        .def("in_", &Pipe::in)  // 'in' is a Python keyword
        .def("out", &Pipe::out)
        .def("all", &Pipe::all)
        .def_ro("cross_section_area", &Pipe::crossSectionArea)
        .def_ro("length", &Pipe::length)
        .def_ro("dispersion_coefficient", &Pipe::dispersion_coefficient)
        .def_ro("cell_distance", &Pipe::cell_distance)
        .def_ro("cell_volume", &Pipe::cell_volume)
        .def("__repr__", [](const Pipe& p) {
            return "<Pipe n_cells=" + std::to_string(p.n_cells) +
                   " length=" + std::to_string(p.length) + ">";
        });

    // ==================== RS_MagneticCaptureProcessChamber ====================

    nb::class_<RS_MagneticCaptureProcessChamber, UnitOperationBase>(
        m, "RS_MagneticCaptureProcessChamber")
        .def(nb::init<const ReactionSystem&, sunindextype, realtype, realtype, realtype,
                      sunindextype, realtype, realtype, realtype, realtype, realtype,
                      std::function<realtype(realtype)>, std::function<realtype(realtype)>,
                      std::function<realtype(realtype)>, std::function<realtype(realtype)>>(),
             "reaction_system"_a, "n_cells_per_phase"_a, "cross_section_area"_a, "length"_a,
             "empty_porosity"_a, "n_disks"_a, "disk_height"_a, "alpha_fluid_particle_volume_ratio"_a,
             "deposition_rate"_a, "MNP_capacity"_a, "magnetic_field_strength"_a,
             "a_eff_function"_a, "capture_function"_a, "flow_rate_function"_a,
             "dispersion_coefficient_function"_a,
             nb::keep_alive<1, 2>())
        .def("set_const_initial_concentration", &RS_MagneticCaptureProcessChamber::setConstInitialConcentration, "concentrations"_a)
        .def("in_", &RS_MagneticCaptureProcessChamber::in)
        .def("out", &RS_MagneticCaptureProcessChamber::out)
        .def("liquid", &RS_MagneticCaptureProcessChamber::liquid)
        .def("slurry", &RS_MagneticCaptureProcessChamber::slurry)
        .def_ro("n_cells_per_phase", &RS_MagneticCaptureProcessChamber::n_cells_per_phase)
        .def_ro("cross_section_area", &RS_MagneticCaptureProcessChamber::crossSectionArea)
        .def_ro("length", &RS_MagneticCaptureProcessChamber::length)
        .def_ro("empty_porosity", &RS_MagneticCaptureProcessChamber::empty_porosity)
        .def("__repr__", [](const RS_MagneticCaptureProcessChamber& pc) {
            return "<RS_MagneticCaptureProcessChamber n_cells=" + std::to_string(pc.n_cells) + ">";
        });

    // ==================== Process ====================

    nb::class_<Process>(m, "Process")
        .def(nb::init<const ComponentSystem&,
                      const std::vector<std::shared_ptr<UnitOperationBase>>&,
                      realtype>(),
             "component_system"_a, "unit_operations"_a,
             "t_end"_a = std::numeric_limits<realtype>::infinity(),
             nb::keep_alive<1, 2>(), nb::keep_alive<1, 3>())
        .def("add_connection", &Process::addConnection,
             "from_mapper"_a, "to_mapper"_a, "flow_rate_function"_a)
        .def("get_unit_operations", &Process::getUnitOperations, nb::rv_policy::reference)
        .def("__repr__", [](const Process& p) {
            return "<Process with " + std::to_string(p.getUnitOperations().size()) + " unit operations>";
        });

    // ==================== SnapshotObserver ====================

    nb::class_<SnapshotObserver>(m, "SnapshotObserver")
        .def(nb::init<realtype, const ArrayMapper&, bool>(),
             "time"_a, "mapper"_a, "compute_errors"_a = false)
        .def_rw("t_desired", &SnapshotObserver::t_desired)
        .def_rw("t_measured", &SnapshotObserver::t_measured)
        .def_ro("mapper", &SnapshotObserver::mapper)
        .def_rw("snapshot", &SnapshotObserver::snapshot)
        .def_rw("compute_errors", &SnapshotObserver::compute_errors)
        .def_rw("error", &SnapshotObserver::error)
        .def("__repr__", [](const SnapshotObserver& o) {
            return "<SnapshotObserver t=" + std::to_string(o.t_desired) + ">";
        });

    // ==================== TimeSeriesObserver ====================

    nb::class_<TimeSeriesObserver>(m, "TimeSeriesObserver")
        .def(nb::init<realtype, realtype, std::size_t, const ArrayMapper&, bool>(),
             "t_start"_a, "t_end"_a, "n_snapshots"_a, "mapper"_a, "compute_errors"_a = false)
        .def(nb::init<realtype, realtype, std::size_t, const ArrayMapper&, Solver&, bool, bool>(),
             "t_start"_a, "t_end"_a, "n_snapshots"_a, "mapper"_a, "solver"_a,
             "add_to_solver"_a = true, "compute_errors"_a = false,
             nb::keep_alive<1, 5>())
        .def("add_all_snapshots_to_solver", &TimeSeriesObserver::add_all_snapshots_to_solver, "solver"_a)
        .def("get_all_snapshots_as_vector", &TimeSeriesObserver::get_all_snapshots_as_vector)
        .def("save_to_npz", &TimeSeriesObserver::save_to_npz, "zipname"_a, "fname"_a, "mode"_a)
        .def("save_to_npy", &TimeSeriesObserver::save_to_npy, "filename"_a)
        .def("__repr__", [](const TimeSeriesObserver&) {
            return "<TimeSeriesObserver>";
        });

    // ==================== Solver ====================

    nb::class_<Solver>(m, "Solver")
        .def(nb::init<const Process&, SolverType>(),
             "process"_a, "solver_type"_a = SolverType::BDF,
             nb::keep_alive<1, 2>())
        .def("solve", &Solver::solve,
             "t_stop"_a = std::numeric_limits<realtype>::infinity(),
             "timeout_seconds"_a = std::numeric_limits<realtype>::infinity())
        .def("get_y_size", &Solver::getYSize)
        .def("get_y", &Solver::getY)
        .def("get_t", &Solver::getT)
        .def("get_internal_time_stamps", &Solver::getInternalTimeStamps)
        .def("log_statistics", &Solver::logStatistics)
        .def("add", &Solver::add, "observer"_a, nb::keep_alive<1, 2>())
        .def("get_solve_time", &Solver::getSolveTime)
        .def("__repr__", [](const Solver& s) {
            return "<Solver t=" + std::to_string(s.getT()) +
                   " y_size=" + std::to_string(s.getYSize()) + ">";
        });

    // ==================== Convenience: Eigen Array types ====================

    // Allow Python to work with RowVector (1xN array)
    // nanobind's eigen support handles this automatically for basic cases

    // ==================== Utility Functions ====================

    m.def("convert_to_pH", &convert_to_pH,
          "observer"_a, "cell_index"_a, "cell_volume"_a, "component_system"_a, "activity_model"_a,
          "Convert TimeSeriesObserver data to pH values for a specific cell");
}
