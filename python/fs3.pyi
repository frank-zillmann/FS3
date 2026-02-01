"""
FS3 - Fast and Flexible Framework for Simple Simulations of Separation-Processes

Python bindings for the FS3 C++ library.
"""

from __future__ import annotations
from typing import Callable, List, Optional, Tuple
import numpy as np
from numpy.typing import NDArray

# Type alias for realtype (double)
realtype = float

# ==================== Enums ====================

class ComponentType:
    """Component type classification."""

    NonMagneticComponent: ComponentType
    MagneticNanoParticle: ComponentType
    MagneticNanoParticleGroup: ComponentType

class Phase:
    """Phase enumeration for multi-phase systems."""

    Liquid: Phase
    SolidOrSlurry: Phase

class UnitOperationType:
    """Unit operation types in the process."""

    Unknown: UnitOperationType
    Inlet: UnitOperationType
    Outlet: UnitOperationType
    Volume: UnitOperationType
    Pipe: UnitOperationType
    RS_MagneticCaptureProcessChamber: UnitOperationType

class SolverType:
    """Available solver types for ODE integration."""

    BDF: SolverType
    ADAMS: SolverType
    ERK: SolverType
    ARK: SolverType

# ==================== Components ====================

class Component:
    """
    Chemical or particulate component metadata.

    Holds name, charge, type and optional physical parameters.
    Provides a fluent builder-style API for concise setup.
    """

    name: str
    charge: int
    type: ComponentType
    molar_mass: float
    radius: Optional[float]
    density: Optional[float]
    magnetic_saturation: Optional[float]
    truesdell_jones_alpha: Optional[float]
    truesdell_jones_beta: Optional[float]

    def __init__(
        self,
        name: str,
        *,
        molar_mass_kg_per_mol: Optional[float] = None,
        charge: Optional[int] = None,
        type: Optional[ComponentType] = None,
        truesdell_jones_alpha: Optional[float] = None,
        truesdell_jones_beta: Optional[float] = None,
        radius_m: Optional[float] = None,
        density_kg_per_m3: Optional[float] = None,
        magnetic_saturation_A_per_m: Optional[float] = None,
    ) -> None: ...

class ComponentSystem:
    """
    Registry of components and medium properties.

    Owns the ordered list of components; caches helper vectors.
    """

    n_components: int
    components: List[Component]
    medium_name: str
    temperature: float
    density: float
    dynamic_viscosity: float
    relative_permittivity: float

    def __init__(self, components: List[Component]) -> None: ...
    def get_idx(self, name: str) -> int: ...
    def get_molar_masses(self) -> NDArray[np.float64]: ...
    def get_inv_molar_masses(self) -> NDArray[np.float64]: ...
    def get_charges(self) -> NDArray[np.float64]: ...
    def __getitem__(self, name_or_index: str | int) -> Component: ...
    def __len__(self) -> int: ...

# ==================== Activity Models ====================

class ActivityModelBase:
    """Base class for activity coefficient models."""

    pass

class NoActivityModel(ActivityModelBase):
    """Identity activity model (activities = concentrations)."""

    def __init__(self, component_system: ComponentSystem) -> None: ...

class TruesdellJonesActivityModel(ActivityModelBase):
    """Truesdell-Jones activity model for electrolytes."""

    def __init__(self, component_system: ComponentSystem) -> None: ...

# ==================== Reactions ====================

class Reaction:
    """Vectorized reaction wrapper."""

    pass

class ReactionSystem:
    """Collection of reactions with an activity model."""

    reactions: List[Reaction]

    def __init__(
        self, component_system: ComponentSystem, activity_model: ActivityModelBase
    ) -> None: ...
    def add(self, reaction: Reaction) -> None: ...
    def __len__(self) -> int: ...

def parse_mass_action_law(
    component_system: ComponentSystem, reaction_str: str
) -> Tuple[List[int], List[int]]:
    """Parse a reaction string and return (reactant_indices, product_indices)."""
    ...

def mass_action_law(
    component_system: ComponentSystem,
    reaction_str: str,
    k_forward: float,
    k_backward: float,
    error_component_idx: int = -1,
) -> Reaction:
    """Create a mass action law reaction from a string like 'H2O <-> H+ + OH-'."""
    ...

def mass_action_law_damped(
    component_system: ComponentSystem,
    reaction_str: str,
    k_forward: float,
    k_backward: float,
    tau_reaction: float,
    error_component_idx: int = -1,
) -> Reaction:
    """Create a damped mass action law reaction with rate smoothing."""
    ...

def mass_action_law_inverse_rate_prediction(
    component_system: ComponentSystem,
    reaction_str: str,
    K_eq: float,
    tau_reaction: float,
    error_component_idx: int = -1,
) -> Reaction:
    """Create an inverse rate prediction mass action law reaction."""
    ...

# ==================== ArrayMapper ====================

class ArrayMapper:
    """Lightweight mapper for strided Array views."""

    rows: int
    cols: int
    start_idx: int

    def __init__(self) -> None: ...

# ==================== Unit Operations ====================

class UnitOperationBase:
    """Abstract base for all unit operations."""

    n_cells: int
    y: NDArray[np.float64]

    def get_type(self) -> UnitOperationType: ...
    def n_components(self) -> int: ...
    def y_size(self) -> int: ...

class Inlet(UnitOperationBase):
    """Boundary inlet unit operation."""

    def __init__(
        self,
        reaction_system: ReactionSystem,
        solution_function: Callable[[float], NDArray[np.float64]],
    ) -> None: ...
    def out(self) -> ArrayMapper: ...
    def get_solution(self, time: float) -> NDArray[np.float64]: ...

class Outlet(UnitOperationBase):
    """Boundary outlet unit operation."""

    def __init__(self, reaction_system: ReactionSystem) -> None: ...
    def in_(self) -> ArrayMapper: ...

class Volume(UnitOperationBase):
    """Single-cell volume unit operation."""

    cell_volume: float

    def __init__(
        self, reaction_system: ReactionSystem, volume: float = float("nan")
    ) -> None: ...
    def set_const_initial_concentration(
        self, concentrations: NDArray[np.float64]
    ) -> None: ...
    def in_(self) -> ArrayMapper: ...
    def out(self) -> ArrayMapper: ...
    def all(self) -> ArrayMapper: ...

class Pipe(UnitOperationBase):
    """Uniform single-phase pipe unit operation."""

    cross_section_area: float
    length: float
    dispersion_coefficient: float
    cell_distance: float
    cell_volume: float

    def __init__(
        self,
        reaction_system: ReactionSystem,
        n_cells: int,
        cross_section_area: float,
        length: float,
        flow_rate_function: Callable[[float], float],
        dispersion_coefficient: float,
    ) -> None: ...
    def set_const_initial_concentration(
        self, concentrations: NDArray[np.float64]
    ) -> None: ...
    def in_(self) -> ArrayMapper: ...
    def out(self) -> ArrayMapper: ...
    def all(self) -> ArrayMapper: ...

class RS_MagneticCaptureProcessChamber(UnitOperationBase):
    """Two-phase rotor-stator magnetic capture chamber."""

    n_cells_per_phase: int
    cross_section_area: float
    length: float
    empty_porosity: float

    def __init__(
        self,
        reaction_system: ReactionSystem,
        n_cells_per_phase: int,
        cross_section_area: float,
        length: float,
        empty_porosity: float,
        n_disks: int,
        disk_height: float,
        alpha_fluid_particle_volume_ratio: float,
        deposition_rate: float,
        MNP_capacity: float,
        magnetic_field_strength: float,
        a_eff_function: Callable[[float], float],
        capture_function: Callable[[float], float],
        flow_rate_function: Callable[[float], float],
        dispersion_coefficient_function: Callable[[float], float],
    ) -> None: ...
    def set_const_initial_concentration(
        self, concentrations: NDArray[np.float64]
    ) -> None: ...
    def in_(self) -> ArrayMapper: ...
    def out(self) -> ArrayMapper: ...
    def liquid(self) -> ArrayMapper: ...
    def slurry(self) -> ArrayMapper: ...

# ==================== Process ====================

class Process:
    """Aggregates unit operations and flow connections."""

    def __init__(
        self,
        component_system: ComponentSystem,
        unit_operations: List[UnitOperationBase],
    ) -> None: ...
    def add_connection(
        self,
        from_mapper: ArrayMapper,
        to_mapper: ArrayMapper,
        flow_rate_function: Callable[[float], float],
    ) -> None: ...
    def get_unit_operations(self) -> List[UnitOperationBase]: ...

# ==================== Observers ====================

class SnapshotObserver:
    """Observer that captures a snapshot at a specific time."""

    t_desired: float
    t_measured: float
    mapper: ArrayMapper
    snapshot: NDArray[np.float64]
    compute_errors: bool
    error: float

    def __init__(
        self, time: float, mapper: ArrayMapper, compute_errors: bool = False
    ) -> None: ...

class TimeSeriesObserver:
    """Observer that captures multiple snapshots over time."""

    def __init__(
        self,
        t_start: float,
        t_end: float,
        n_snapshots: int,
        mapper: ArrayMapper,
        compute_errors: bool = False,
    ) -> None: ...
    def add_all_snapshots_to_solver(self, solver: Solver) -> None: ...
    def get_all_snapshots_as_vector(
        self,
    ) -> Tuple[List[float], Tuple[int, int, int]]: ...
    def save_to_npz(self, zipname: str, fname: str, mode: str) -> None: ...
    def save_to_npy(self, filename: str) -> None: ...

# ==================== Solver ====================

class Solver:
    """SUNDIALS-based ODE solver wrapper."""

    def __init__(
        self, process: Process, solver_type: SolverType = SolverType.BDF
    ) -> None: ...
    def solve(
        self, t_stop: float = float("inf"), timeout_seconds: float = float("inf")
    ) -> None: ...
    def get_y_size(self) -> int: ...
    def get_y(self) -> List[float]: ...
    def get_t(self) -> float: ...
    def get_internal_time_stamps(self) -> List[float]: ...
    def log_statistics(self) -> None: ...
    def add(self, observer: SnapshotObserver) -> None: ...
    def get_solve_time(self) -> float: ...
