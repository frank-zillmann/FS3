#ifndef COMPONENT_HPP
#define COMPONENT_HPP

#include <optional>
#include <string>

#include "sundials/sundials_types.h"

/**
 * @brief Classification of a component's physical nature.
 */
enum class Type { NonMagneticComponent, MagneticNanoParticle, MagneticNanoParticleGroup };

/**
 * @brief Represents a single chemical or physical species in the system.
 *
 * Stores name, charge, type, and optional physical parameters (molar mass,
 * activity coefficients, particle geometry, magnetic saturation).
 * Supports a fluent builder interface for convenient construction.
 */
class Component {
   public:
    std::string name;
    int charge = 0;                          ///< Elementary charge (e.g. +1 for H⁺, -1 for OH⁻)
    Type type = Type::NonMagneticComponent;  ///< Component classification
    realtype molarMass = 1;                  ///< Molar mass [g/mol]
    std::optional<realtype> truesdell_jones_alpha = std::nullopt;  ///< Truesdell-Jones ion-size parameter [m]
    std::optional<realtype> truesdell_jones_beta = std::nullopt;  ///< Truesdell-Jones salting-out coefficient [m³/mol]
    std::optional<realtype> radius = std::nullopt;                ///< Particle radius [m]
    std::optional<realtype> density = std::nullopt;               ///< Particle density [kg/m³]
    std::optional<realtype> magnetic_saturation = std::nullopt;   ///< Magnetic saturation [A/m]

    /// Constructs a Component with the given name; all other fields use defaults.
    explicit Component(const std::string& name) : name(name) {}

    Component& setCharge(int integer_charge);
    Component& setType(Type component_type);

    /// @throws std::invalid_argument if molarMass_kg_per_mol <= 0
    Component& setMolarMass(realtype molarMass_kg_per_mol);

    /// Sets both Truesdell-Jones parameters; warns if the result degenerates to Debye-Hückel.
    /// @throws std::invalid_argument if either parameter is negative
    Component& setTruesdellJonesParameters(realtype alpha_meter = 0.0, realtype beta_cubicmeter_per_mol = 0.0);

    /// @throws std::invalid_argument if radius_m <= 0
    Component& setRadius(realtype radius_m);

    /// @throws std::invalid_argument if density_kg_per_m3 <= 0
    Component& setDensity(realtype density_kg_per_m3);

    /// @throws std::invalid_argument if magnetic_saturation_A_per_m <= 0
    Component& setMagneticSaturation(realtype magnetic_saturation_A_per_m);
};

#endif  // COMPONENT_HPP