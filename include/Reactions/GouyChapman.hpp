#ifndef GOUY_CHAPMAN_HPP
#define GOUY_CHAPMAN_HPP

#include <cstddef>
#include <string>

#include "Components/ComponentSystem.hpp"
#include "EigenDataTypes.hpp"

class GouyChapmanModel {
   public:
    GouyChapmanModel(const ComponentSystem& componentSystem,
                     const std::string& mnp_oh2_plus_name,
                     const std::string& mnp_oh_name,
                     const std::string& mnp_o_minus_name,
                     const std::string& h_plus_name,
                     realtype surface_site_density_mol_per_m2,
                     realtype ion_concentration_floor = 1e-5,
                     realtype surface_group_epsilon = 1e-12);

    void operator()(realtype t, ArrayMap& concentrations, ArrayMap& activities, ArrayMap& dc_dt) const;

   private:
    std::size_t idx_mnp_oh2_plus;
    std::size_t idx_mnp_oh;
    std::size_t idx_mnp_o_minus;
    std::size_t idx_h_plus;
    RowVector charges_nonzero_mask;
    realtype surface_charge_factor;
    realtype normalization_factor;
    realtype ion_concentration_floor;
    realtype surface_group_epsilon;
};

#endif  // GOUY_CHAPMAN_HPP