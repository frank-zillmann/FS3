// Langmuir binding reaction helper
#ifndef LANGMUIR_HPP
#define LANGMUIR_HPP

#include <cstddef>
#include <functional>

#include "EigenDataTypes.hpp"
#include "Reactions/Reaction.hpp"

// K_B callback returns L/g = m^3/kg; converted internally to m^3/mol.
Reaction langmuir_binding_reaction(
	std::size_t idx_h_plus,
	std::size_t idx_mnp1,
	std::size_t idx_mnp10,
	std::size_t idx_mnp_higg,
	std::size_t idx_higg,
	realtype molar_mass_higg,
	std::function<realtype(realtype)> k_b_from_pH,
	std::function<realtype(realtype)> q_max_from_pH,
	realtype tau_reaction);

#endif  // LANGMUIR_HPP
