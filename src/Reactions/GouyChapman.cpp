#include "Reactions/GouyChapman.hpp"

#include <sstream>

#include "Logger.hpp"
#include "PhysicalConstants.hpp"

GouyChapmanModel::GouyChapmanModel(const ComponentSystem& componentSystem,
                                   const std::string& mnp_oh2_plus_name,
                                   const std::string& mnp_oh_name,
                                   const std::string& mnp_o_minus_name,
                                   const std::string& h_plus_name,
                                   realtype surface_site_density_mol_per_m2,
                                   realtype ion_concentration_floor,
                                   realtype surface_group_epsilon)
    : idx_mnp_oh2_plus(componentSystem.getIdx(mnp_oh2_plus_name)),
      idx_mnp_oh(componentSystem.getIdx(mnp_oh_name)),
      idx_mnp_o_minus(componentSystem.getIdx(mnp_o_minus_name)),
      idx_h_plus(componentSystem.getIdx(h_plus_name)),
      charges_nonzero_mask((componentSystem.getCharges() != 0).cast<realtype>()),
      surface_charge_factor(constants::elementary_charge * constants::avogadro * surface_site_density_mol_per_m2),
      normalization_factor(2 * componentSystem.relative_permittivity * constants::vacuum_permittivity *
                           constants::gas_constant * componentSystem.temperature),
      ion_concentration_floor(ion_concentration_floor),
      surface_group_epsilon(surface_group_epsilon) {}

void GouyChapmanModel::operator()(realtype t, ArrayMap& concentrations, ArrayMap& activities, ArrayMap& dc_dt) const {
    (void)t;
    (void)dc_dt;
    const auto c_mnp_oh2_plus = concentrations.col(idx_mnp_oh2_plus);
    const auto c_mnp_oh = concentrations.col(idx_mnp_oh);
    const auto c_mnp_o_minus = concentrations.col(idx_mnp_o_minus);
    const auto c_h_plus = concentrations.col(idx_h_plus);

    const auto surface_groups_total = c_mnp_oh2_plus + c_mnp_oh + c_mnp_o_minus + surface_group_epsilon;
    const auto surface_charge = surface_charge_factor * (c_mnp_oh2_plus - c_mnp_o_minus) / surface_groups_total;

    const auto ion_concentration = (concentrations.rowwise() * charges_nonzero_mask).rowwise().sum() +
                                   ion_concentration_floor;

    const auto normalized_surface_charge = 0.5 * (surface_charge / ((normalization_factor * ion_concentration).sqrt()));

    const auto f = (-normalized_surface_charge + (normalized_surface_charge.square() + 1).sqrt()).square();

    const auto c_h_plus_surface = f * c_h_plus;

    activities.col(idx_h_plus) = c_h_plus_surface;

#if LOG_ENABLED
    static int call_count = 0;
    if (call_count < LOG_FIRST_N_CALLS || call_count % LOG_EVERY_N_CALLS == 0) {
        std::ostringstream oss;
        oss << "Call count: " << call_count << " at t= " << t << "\n";
        oss << "H+ concentration: " << c_h_plus.transpose() << "\n";
        oss << "MNP-OH2+ concentration: " << c_mnp_oh2_plus.transpose() << "\n";
        oss << "MNP-OH concentration : " << c_mnp_oh.transpose() << "\n";
        oss << "MNP-O- concentration: " << c_mnp_o_minus.transpose() << "\n";
        oss << "Surface charge : " << surface_charge.transpose() << "\n";
        oss << "Ion concentration: " << ion_concentration.transpose() << "\n";
        oss << "Normalized surface charge : " << normalized_surface_charge.transpose() << "\n";
        oss << "f: " << f.transpose() << "\n\n";
        LOG("MNP_H+_gouy_chapman.log", oss.str());
    }
    call_count++;
#endif
}
