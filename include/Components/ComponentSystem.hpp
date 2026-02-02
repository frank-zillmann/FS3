#ifndef COMPONENTSYSTEM_HPP
#define COMPONENTSYSTEM_HPP

#include <vector>

#include "Components/Component.hpp"
#include "EigenDataTypes.hpp"
#include "sundials/sundials_types.h"

/**
 * @brief Registry of components and medium properties
 *
 * Owns the ordered list of components; caches helper vectors (molar masses,
 * inverse molar masses, charges). Provides mass⟷molar concentration conversion
 * utilities and name/idx lookup.
 */
class ComponentSystem {
   public:
    const std::vector<Component> components;
    const sunindextype n_components;

    // Medium/Conditions the components are in:
    const realtype temperature;                    // [K]
    const realtype density;                        // [kg/m^3]
    const realtype dynamic_viscosity;              // [Pa*s]
    const realtype relative_permittivity = 78.54;  // [1] (dimensionless)

    ComponentSystem(const std::vector<Component>& components,
                    realtype temperature = 298.15,
                    realtype density = 1000.0,
                    realtype dynamic_viscosity = 0.001,
                    realtype relative_permittivity = 78.54)
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

    std::size_t getIdx(const std::string& name) const;

    const Component& operator()(const std::string& name) const;

    const Component& operator()(sunindextype index) const;

    template <EigenArrayLike InputArrayType, WritableEigenArrayLike OutputArrayType>
    void massConcentrationToMolarConcentration(const InputArrayType& mass_concentrations,
                                               OutputArrayType& molar_concentrations) const;

    // Convenience overload that returns a new Array
    template <EigenArrayLike InputArrayType>
    auto massConcentrationToMolarConcentration(const InputArrayType& mass_concentrations) const -> Array;

    template <EigenArrayLike InputArrayType, WritableEigenArrayLike OutputArrayType>
    void molarConcentrationToMassConcentration(const InputArrayType& molar_concentrations,
                                               OutputArrayType& mass_concentrations) const;

    // Convenience overload that returns a new Array
    template <EigenArrayLike InputArrayType>
    auto molarConcentrationToMassConcentration(const InputArrayType& molar_concentrations) const -> Array;

    const RowVector& getMolarMasses() const { return molar_masses; }
    const RowVector& getCharges() const { return charges; }

   private:
    // Cached vectors for performance
    RowVector molar_masses;
    RowVector charges;
    void initializeHelperVectors();
};

// Include template implementations
#include "Components/ComponentSystem.tpp"

#endif  // COMPONENTSYSTEM_HPP