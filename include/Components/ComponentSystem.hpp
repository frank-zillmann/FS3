#ifndef COMPONENTSYSTEM_HPP
#define COMPONENTSYSTEM_HPP

#include <vector>

#include "Components/Component.hpp"
#include "EigenDataTypes.hpp"
#include "sundials/sundials_types.h"

/**
 * @brief Owns the ordered list of components and ambient medium properties.
 *
 * Provides lookup by name and index, unit-conversion helpers between mass and
 * molar concentrations, and cached helper vectors for performance-critical code.
 */
class ComponentSystem {
   public:
    const std::vector<Component> components;
    const sunindextype n_components;

    // Ambient medium properties:
    const realtype temperature;            ///< [K]
    const realtype density;                ///< [kg/m³]
    const realtype dynamic_viscosity;      ///< [Pa·s]
    const realtype relative_permittivity;  ///< [–] (dimensionless)

    ComponentSystem(const std::vector<Component>& components,
                    realtype temperature = 298.15,
                    realtype density = 1000.0,
                    realtype dynamic_viscosity = 0.001,
                    realtype relative_permittivity = 78.54);

    /// Returns the index of the component with the given name.
    /// @throws std::runtime_error if the name is not found.
    std::size_t getIdx(const std::string& name) const;

    /// Returns the component with the given name.
    /// @throws std::runtime_error if the name is not found.
    const Component& operator()(const std::string& name) const;

    /// Returns the component at the given index.
    /// @throws std::out_of_range if the index is out of bounds.
    const Component& operator()(sunindextype index) const;

    /// Converts mass concentrations [g/m³] to molar concentrations [mol/m³] in-place.
    template <EigenArrayLike InputArrayType, WritableEigenArrayLike OutputArrayType>
    void massConcentrationToMolarConcentration(const InputArrayType& mass_concentrations,
                                               OutputArrayType& molar_concentrations) const;

    /// Convenience overload returning a new Array.
    template <EigenArrayLike InputArrayType>
    auto massConcentrationToMolarConcentration(const InputArrayType& mass_concentrations) const -> Array;

    /// Converts molar concentrations [mol/m³] to mass concentrations [g/m³] in-place.
    template <EigenArrayLike InputArrayType, WritableEigenArrayLike OutputArrayType>
    void molarConcentrationToMassConcentration(const InputArrayType& molar_concentrations,
                                               OutputArrayType& mass_concentrations) const;

    /// Convenience overload returning a new Array.
    template <EigenArrayLike InputArrayType>
    auto molarConcentrationToMassConcentration(const InputArrayType& molar_concentrations) const -> Array;

    const RowVector& getMolarMasses() const { return molar_masses; }
    const RowVector& getCharges() const { return charges; }

   private:
    RowVector molar_masses;  ///< Cached molar masses [g/mol], indexed by component order
    RowVector charges;       ///< Cached charges [–], indexed by component order
    void initializeHelperVectors();
};

// Include template implementations
#include "Components/ComponentSystem.tpp"

#endif  // COMPONENTSYSTEM_HPP