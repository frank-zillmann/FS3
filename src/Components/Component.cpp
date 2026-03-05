#include "Components/Component.hpp"

#include <iostream>
#include <stdexcept>

Component& Component::setCharge(int integer_charge) {
    charge = integer_charge;
    return *this;
}

Component& Component::setType(Type component_type) {
    type = component_type;
    return *this;
}

Component& Component::setMolarMass(realtype molarMass_kg_per_mol) {
    if (molarMass_kg_per_mol <= 0.0) {
        throw std::invalid_argument("Molar mass must be positive!");
    }
    molarMass = molarMass_kg_per_mol;
    return *this;
}

Component& Component::setTruesdellJonesParameters(realtype alpha_meter, realtype beta_cubicmeter_per_mol) {
    truesdell_jones_alpha = alpha_meter;
    truesdell_jones_beta = beta_cubicmeter_per_mol;
    if (alpha_meter < 0.0 || beta_cubicmeter_per_mol < 0.0) {
        throw std::invalid_argument("Truesdell-Jones parameters must be non-negative!");
    } else if (alpha_meter == 0.0 && beta_cubicmeter_per_mol == 0.0) {
        std::cerr << "Warning: Truesdell-Jones parameters are both zero, activity correction behaves according to "
                     "Debye-Hückel model."
                  << std::endl;
    } else if (alpha_meter == 0) {
        std::cerr << "Warning: Truesdell-Jones parameter alpha_nm is zero, activity correction behaves according to "
                     "Extended Debye-Hückel model."
                  << std::endl;
    }
    return *this;
}

Component& Component::setRadius(realtype radius_m) {
    if (radius_m <= 0.0) {
        throw std::invalid_argument("Radius must be positive!");
    }
    radius = radius_m;
    return *this;
}

Component& Component::setDensity(realtype density_kg_per_m3) {
    if (density_kg_per_m3 <= 0.0) {
        throw std::invalid_argument("Density must be positive!");
    }
    density = density_kg_per_m3;
    return *this;
}

Component& Component::setMagneticSaturation(realtype magnetic_saturation_A_per_m) {
    if (magnetic_saturation_A_per_m <= 0.0) {
        throw std::invalid_argument("Magnetic saturation must be positive!");
    }
    magnetic_saturation = magnetic_saturation_A_per_m;
    return *this;
}