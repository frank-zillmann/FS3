#include "Components/ComponentSystem.hpp"

#include <algorithm>
#include <stdexcept>

ComponentSystem::ComponentSystem(const std::vector<Component>& components,
                                 realtype temperature,
                                 realtype density,
                                 realtype dynamic_viscosity,
                                 realtype relative_permittivity)
    : components(components),
      n_components(components.size()),
      molar_masses(components.size()),
      charges(components.size()),
      temperature(temperature),
      density(density),
      dynamic_viscosity(dynamic_viscosity),
      relative_permittivity(relative_permittivity) {
    initializeHelperVectors();
}

std::size_t ComponentSystem::getIdx(const std::string& name) const {
    auto it = std::find_if(components.begin(), components.end(),
                           [&name](const Component& comp) { return comp.name == name; });
    if (it != components.end()) {
        return std::distance(components.begin(), it);
    } else {
        throw std::runtime_error("Component not found: " + name);
    }
}

const Component& ComponentSystem::operator()(const std::string& name) const {
    auto it = std::find_if(components.begin(), components.end(),
                           [&name](const Component& comp) { return comp.name == name; });
    if (it != components.end()) {
        return *it;
    } else {
        throw std::runtime_error("Component not found: " + name);
    }
}

const Component& ComponentSystem::operator()(sunindextype index) const {
    if (index >= 0 && static_cast<std::size_t>(index) < components.size()) {
        return components[index];
    } else {
        throw std::out_of_range("Index out of range: " + std::to_string(index));
    }
}

void ComponentSystem::initializeHelperVectors() {
    for (sunindextype i = 0; i < n_components; ++i) {
        molar_masses(i) = components[i].molarMass;
        charges(i) = components[i].charge;
    }
}
