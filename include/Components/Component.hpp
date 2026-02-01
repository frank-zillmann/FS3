#ifndef COMPONENT_HPP
#define COMPONENT_HPP

#include <cassert>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>

#include "sundials/sundials_types.h"

/**
 * @brief Component type classification
 */
enum class Type { NonMageneticComponent, MagneticNanoParticle, MagneticNanoParticleGroup };

/**
 * @brief Chemical or particulate component metadata
 *
 * Holds name, charge, type and optional physical parameters (molar mass,
 * particle radius/density, magnetic saturation, activity parameters).
 * Provides a fluent builder-style API for concise setup.
 */
class Component {
   public:
    std::string name;
    int charge = 0;                           // Charge of the component, e.g., +1 for H+, -1 for OH-
    Type type = Type::NonMageneticComponent;  // Type of the component, default is non-magnetic
    realtype molarMass = 1;                   // Molar mass in g/mol
    std::optional<realtype> truesdell_jones_alpha = std::nullopt;  // ion specific parameter for the Truesdell-Jones activity model
    std::optional<realtype> truesdell_jones_beta = std::nullopt;  // ion specific parameter for the Truesdell-Jones activity model
    std::optional<realtype> radius = std::nullopt;                // Radius in m of particles
    std::optional<realtype> density = std::nullopt;               // Density in kg/m^3 of particles
    std::optional<realtype> magnetic_saturation = std::nullopt;  // Magnetic saturation of magnetic particles

    // Constructor that acts as builder entry point
    explicit Component(const std::string& name) : name(name) {}

    // Builder methods that return *this for chaining
    Component& setCharge(int integer_charge) {
        charge = integer_charge;
        return *this;
    }
    Component& setType(Type type) {
        type = type;
        return *this;
    }
    Component& setMolarMass(realtype molarMass_kg_per_mol) {
        if (molarMass_kg_per_mol <= 0.0) {
            throw std::invalid_argument("Molar mass must be positive!");
        }
        molarMass = molarMass_kg_per_mol;
        return *this;
    }
    Component& setTruesdellJonesParameters(realtype alpha_meter = 0.0, realtype beta_cubicmeter_per_mol = 0.0) {
        truesdell_jones_alpha = alpha_meter;
        truesdell_jones_beta = beta_cubicmeter_per_mol;
        if (alpha_meter < 0.0 || beta_cubicmeter_per_mol < 0.0) {
            throw std::invalid_argument("Truesdell-Jones parameters must be non-negative!");
        } else if (alpha_meter == 0.0 && beta_cubicmeter_per_mol == 0.0) {
            std::cerr << "Warning: Truesdell-Jones parameters are both zero, activity correction behaves according to "
                         "Debye-Hückel model."
                      << std::endl;
        } else if (alpha_meter == 0) {
            std::cerr
                << "Warning: Truesdell-Jones parameter alpha_nm is zero, activity correction behaves according to "
                   "Extended Debye-Hückel model."
                << std::endl;
        }
        return *this;
    }
    Component& setRadius(realtype radius_m) {
        radius = radius_m;
        return *this;
    }
    Component& setDensity(realtype density_kg_per_m3) {
        density = density_kg_per_m3;
        return *this;
    }
    Component& setMagneticSaturation(realtype magnetic_saturation_A_per_m) {
        magnetic_saturation = magnetic_saturation_A_per_m;
        return *this;
    }
};

#endif  // COMPONENT_HPP