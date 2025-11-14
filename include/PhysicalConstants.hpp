#ifndef PHYSICAL_CONSTANTS_HPP
#define PHYSICAL_CONSTANTS_HPP

namespace constants {

// Fundamental constants (SI units)
constexpr realtype pi = 3.14159265358979323846;             // π
constexpr realtype avogadro = 6.02214076e23;                // mol⁻¹
constexpr realtype boltzmann = 1.380649e-23;                // J/K
constexpr realtype planck = 6.62607015e-34;                 // J·s
constexpr realtype elementary_charge = 1.602176634e-19;     // C
constexpr realtype gas_constant = 8.314462618;              // J/(mol·K)
constexpr realtype vacuum_permittivity = 8.8541878128e-12;  // F/m (ε₀)
constexpr realtype vacuum_permeability = 1.25663706212e-6;  // H/m (μ₀)

}  // namespace constants

#endif  // PHYSICAL_CONSTANTS_HPP