#include "Reactions/Langmuir.hpp"

Reaction langmuir_binding_reaction(
	std::size_t idx_h_plus,
	std::size_t idx_mnp1,
	std::size_t idx_mnp10,
	std::size_t idx_mnp_higg,
	std::size_t idx_higg,
	realtype molar_mass_higg,
	std::function<realtype(realtype)> k_b_from_pH,
	std::function<realtype(realtype)> q_max_from_pH,
	realtype tau_reaction) {
	return Reaction(
		[idx_h_plus, idx_mnp1, idx_mnp10, idx_mnp_higg, idx_higg, molar_mass_higg, k_b_from_pH, q_max_from_pH,
		 tau_reaction](realtype, const ConstArrayMap& concentrations, const ConstArrayMap& activities, ArrayMap& dc_dt) {
			const auto pH = -(concentrations.col(idx_h_plus)).log10() + 3.0;
			const auto c_mnp = activities.middleCols(idx_mnp1, idx_mnp10 - idx_mnp1 + 1).rowwise().sum();
			const auto c_higg = activities.col(idx_higg);
			const auto c_mnp_higg = activities.col(idx_mnp_higg);

			auto dc_dt_mnp_higg = dc_dt.col(idx_mnp_higg);
			auto dc_dt_higg = dc_dt.col(idx_higg);

			const auto k_b_raw = pH.unaryExpr([&](realtype p) { return k_b_from_pH(p); });
			const auto k_b = k_b_raw * molar_mass_higg;
			const auto q_max = pH.unaryExpr([&](realtype p) { return q_max_from_pH(p); });

			const auto p = c_mnp_higg - c_higg - c_mnp * q_max - molar_mass_higg / k_b;
			const auto q_tilde = c_higg * (c_mnp * q_max - c_mnp_higg) - (molar_mass_higg * c_mnp_higg / k_b);
			const auto radicant = (p.square() - 4.0 * q_tilde).sqrt();
			const auto delta_xi = (-p - radicant) / 2.0;
			const auto net_rate = delta_xi / tau_reaction;

			dc_dt_mnp_higg += net_rate;
			dc_dt_higg -= net_rate;
		});
}
