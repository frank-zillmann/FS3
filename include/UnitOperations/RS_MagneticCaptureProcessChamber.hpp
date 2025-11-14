#ifndef RS_MAGNETIC_CAPTURE_PROCESS_CHAMBER_HPP
#define RS_MAGNETIC_CAPTURE_PROCESS_CHAMBER_HPP

#include <sundials/sundials_types.h>

#include <functional>

#include "EigenDataTypes.hpp"
#include "UnitOperationBase.hpp"

/**
 * @brief Two-phase rotor–stator magnetic capture chamber
 *
 * Discretizes the process chamber with liquid and slurry phases per axial cell.
 * Implements magnetic capture/uncapture of MNPs, inter-phase mass transfer,
 * dynamic volumes (porosity changes), and transport in the liquid phase.
 */
class RS_MagneticCaptureProcessChamber : public UnitOperationBase {
   public:
    RS_MagneticCaptureProcessChamber(const ReactionSystem& reactionSystem,
                                     sunindextype n_cells_per_phase,
                                     realtype crossSectionArea,
                                     realtype length,
                                     realtype empty_porosity,
                                     sunindextype n_disks,
                                     realtype disk_height,
                                     realtype alpha_fluid_particle_volume_ratio,
                                     realtype deposition_rate,
                                     realtype MNP_capacity,
                                     realtype magnetic_field_strength,
                                     std::function<realtype(realtype)> a_eff_function,
                                     std::function<realtype(realtype)> capture_function,
                                     std::function<realtype(realtype)> flowRateFunction,
                                     std::function<realtype(realtype)> dispersion_coefficient_function)
        : UnitOperationBase(reactionSystem, n_cells_per_phase * 2),
          n_cells_per_phase(n_cells_per_phase),
          crossSectionArea(crossSectionArea),
          length(length),
          empty_porosity(empty_porosity),
          n_disks(n_disks),
          disk_height(disk_height),
          alpha_fluid_particle_volume_ratio(alpha_fluid_particle_volume_ratio),
          deposition_rate(deposition_rate),
          MNP_capacity(MNP_capacity),
          magnetic_field_strength(magnetic_field_strength),
          a_eff_function(a_eff_function),
          capture_function(capture_function),
          flowRateFunction(flowRateFunction),
          dispersion_coefficient_function(dispersion_coefficient_function),
          cell_volume_l_and_sl((empty_porosity * crossSectionArea * length) / n_cells_per_phase),
          cell_distance(length / n_cells_per_phase),
          c_l(Array(n_cells_per_phase, n_components())),
          activities_l(Array(n_cells_per_phase, n_components())),
          dc_dt_l(Array(n_cells_per_phase, n_components())),
          c_sl(Array(n_cells_per_phase, n_components())),
          activities_sl(Array(n_cells_per_phase, n_components())),
          dc_dt_sl(Array(n_cells_per_phase, n_components())),
          dm_dt_capture_MNP(Array(n_cells_per_phase, MNP_range.second - MNP_range.first + 1)) {
        init();
    }

    virtual ~RS_MagneticCaptureProcessChamber() = default;

    UnitOperationType getType() const override { return UnitOperationType::RS_MagneticCaptureProcessChamber; }

    void setConstInitialConcentration(const RowVector& concentrations) {
        assert(concentrations.size() == n_components() &&
               "Initial concentration vector size does not match number of components in the component system.");
        y = Array::Zero(n_cells, n_components());
        for (sunindextype i = 0; i < n_cells; i += 2) {
            y.row(i) = cell_volume_l_and_sl * concentrations;  // += 2 to only fill the liquid phase
        }
    };

    const ColVector& get_V_l(realtype t, const realtype* y) const override;

    const ArrayMapper in() const { return ArrayMapper(this, 1, n_components(), 0, 2); }

    const ArrayMapper out() const {
        return ArrayMapper(this, 1, n_components(), 2 * n_components() * (n_cells_per_phase - 1), 2);
    }

    const ArrayMapper liquid() const { return ArrayMapper(this, n_cells_per_phase, n_components(), 0, 2); };

    const ArrayMapper slurry() const { return ArrayMapper(this, n_cells_per_phase, n_components(), n_components(), 2); }

    void rhs_capture(realtype t,
                     const ConstArrayMap& m_l,
                     const ConstArrayMap& m_sl,
                     ArrayMap& dm_dt_l,
                     ArrayMap& dm_dt_sl,
                     realtype flowRate);

    sunindextype idx(sunindextype finiteVolumeIdx, sunindextype componentIdx, Phase phase) const override;

    realtype errorFunction(realtype t, const realtype* y) const override;

    // Override the rhs method to implement the magnetic capture process
    void rhs(realtype t,
             const realtype* y,
             realtype* dy_dt,
             const Solver& solver,
             bool enable_reactions = true,
             bool enable_convection = true,
             bool enable_dispersion = true,
             bool enable_otherPhysics = true) override;

    const sunindextype n_cells_per_phase;  // number liquid cells and slurry cells -> total cells = n_cells_per_phase * 2
    const realtype crossSectionArea;       // Cross-sectional area of the chamber
    const realtype length;                 // Height of the discs in the chamber
    const realtype empty_porosity;         // Porosity of the chamber when empty
    const sunindextype n_disks;            // Number of disks in the chamber
    const realtype disk_height;            // Height of the discs in the chamber
    const realtype alpha_fluid_particle_volume_ratio;  // Ratio of fluid to particle volume in the slurry
    const realtype deposition_rate;                    // Empirical parameter in G to account for deposition rate
    const realtype MNP_capacity;  // Maximum mass concentration of magnetic particles that can be captured
    const realtype magnetic_field_strength;  // In the uncapture case, half of the MNPs in slurry phase decay back to liquid phase every T_halfLife seconds
    const realtype dispersion_coefficient_sl_to_l_NMC = 0.0;  // Dispersion coefficient for NMC between slurry and liquid phases

    const std::function<realtype(realtype)> a_eff_function;  // input is u_m/u_0
    const std::function<realtype(realtype)> capture_function;
    const std::function<realtype(realtype)> flowRateFunction;  // Function to calculate flow rate based on time
    const std::function<realtype(realtype)> dispersion_coefficient_function;  // Function to calculate dispersion coefficient based on time

    // derived helper parameters
    const realtype cell_volume_l_and_sl;  // Volume of each cell (liquid or slurry)
    const realtype cell_distance;         // Distance between the centers of two adjacent cells

   private:
    void init();

    std::pair<sunindextype, sunindextype> NMC_range;
    std::pair<sunindextype, sunindextype> MNP_range;
    std::pair<sunindextype, sunindextype> MNPG_range;
    RowVector radii_MNP;
    RowVector densities_MNP;
    RowVector magnetic_saturations_MNP;

    mutable ColVector V_l;

    mutable Array c_l;
    mutable Array activities_l;
    mutable Array dc_dt_l;
    mutable Array c_sl;
    mutable Array activities_sl;
    mutable Array dc_dt_sl;

    mutable Array dm_dt_capture_MNP;

#if LOG_ENABLED
    int call_count = 0;
#endif
};

#endif  // RS_MAGNETIC_CAPTURE_PROCESS_CHAMBER_HPP