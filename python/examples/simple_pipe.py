"""
Example: Simple pipe simulation using FS3 Python bindings.

This example demonstrates how to:
1. Define components
2. Create a reaction system
3. Set up unit operations (inlet, pipe, outlet)
4. Create a process with connections
5. Run a simulation using the Solver
"""

import numpy as np
import fs3


def main():
    # ==================== Components ====================

    components = [
        fs3.Component("H2O"),
        fs3.Component(
            "H+", charge=1, truesdell_jones_alpha=4.78e-10, truesdell_jones_beta=0.24e-3
        ),
        fs3.Component(
            "OH-",
            charge=-1,
            truesdell_jones_alpha=10.65e-10,
            truesdell_jones_beta=0.21e-3,
        ),
        fs3.Component("Tracer"),
    ]

    component_system = fs3.ComponentSystem(components)
    print(f"Created {component_system}")

    # ==================== Activity Model & Reactions ====================

    activity_model = fs3.NoActivityModel(component_system)
    reaction_system = fs3.ReactionSystem(component_system, activity_model)

    # Add water dissociation reaction: H2O <=> H+ + OH-
    # k_forward and k_backward chosen to give K_eq = 1e-14
    water_reaction = fs3.mass_action_law(
        component_system,
        "H2O <=> H+ + OH-",
        k_forward=1e-3,
        k_backward=1e11,
        error_component_idx=1,  # Monitor H+ for pH calculation
    )
    reaction_system.add(water_reaction)
    print(f"Created {reaction_system}")

    # ==================== Unit Operations ====================

    # Inlet: provides boundary conditions
    def inlet_concentration(t: float) -> np.ndarray:
        """Return inlet concentrations as function of time."""
        # [H2O, H+, OH-, Tracer] in mol/m³ (= mM)
        tracer = 1.0 if t < 10.0 else 0.0  # Pulse injection
        return np.array([55500.0, 1e-4, 1e-4, tracer])

    inlet = fs3.Inlet(reaction_system, inlet_concentration)

    # Pipe: 10 cells, 1 cm² cross-section, 10 cm length
    flow_rate = lambda t: 1e-6  # 1 mL/s = 1e-6 m³/s
    pipe = fs3.Pipe(
        reaction_system,
        n_cells=10,
        cross_section_area=1e-4,  # 1 cm² = 1e-4 m²
        length=0.1,  # 10 cm = 0.1 m
        flow_rate_function=flow_rate,
        dispersion_coefficient=1e-8,
    )

    # Set initial conditions (same as inlet)
    pipe.set_const_initial_concentration(np.array([55500.0, 1e-4, 1e-4, 0.0]))

    # Outlet: collects outflow
    outlet = fs3.Outlet(reaction_system)

    print(f"Created {inlet}")
    print(f"Created {pipe}")
    print(f"Created {outlet}")

    # ==================== Process ====================

    process = fs3.Process(component_system, [inlet, pipe, outlet])

    # Connect: inlet -> pipe -> outlet
    process.add_connection(inlet.out(), pipe.in_(), flow_rate)
    process.add_connection(pipe.out(), outlet.in_(), flow_rate)

    print(f"Created {process}")

    # ==================== Solver ====================

    solver = fs3.Solver(process, fs3.SolverType.BDF)

    print(f"Created {solver}")

    t_end = 30.0  # Simulate for 30 seconds

    # Add observers to capture results at pipe outlet
    observer = fs3.TimeSeriesObserver(
        t_start=0.0,
        t_end=t_end,
        n_snapshots=100,
        mapper=pipe.out(),
        compute_errors=False,
    )
    observer.add_all_snapshots_to_solver(solver)

    print(f"Added observer: {observer}")

    # ==================== Simulation ====================
    print("Running simulation...")

    solver.solve(t_stop=t_end, timeout_seconds=60.0)

    print(f"Simulation completed!")
    print(f"  Final time: {solver.get_t():.3f} s")
    print(f"  Solve time: {solver.get_solve_time():.3f} s")
    print(f"  Number of internal steps: {len(solver.get_internal_time_stamps())}")

    # Get results
    data, shape = observer.get_all_snapshots_as_vector()
    results = np.array(data).reshape(shape)
    print(f"  Results shape: {results.shape}  # (n_snapshots, n_cells, n_components)")

    # Save results
    observer.save_to_npy("pipe_simulation_results.npy")
    print("  Results saved to pipe_simulation_results.npy")

    import matplotlib.pyplot as plt

    # Plot tracer concentration at pipe outlet over time
    time_points = np.linspace(0, t_end, results.shape[0])
    tracer_conc_outlet = results[:, -1, 3]  # Last cell, Tracer component
    plt.plot(time_points, tracer_conc_outlet)
    plt.xlabel("Time (s)")
    plt.ylabel("Tracer Concentration (mol/m³)")
    plt.title("Tracer Concentration at Pipe Outlet Over Time")
    plt.grid()
    plt.show()


if __name__ == "__main__":
    main()
