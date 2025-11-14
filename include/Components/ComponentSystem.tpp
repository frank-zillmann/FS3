#ifndef COMPONENTSYSTEM_TPP
#define COMPONENTSYSTEM_TPP

#include <cassert>

// Template implementations for ComponentSystem

template <EigenArrayLike InputArrayType, WritableEigenArrayLike OutputArrayType>
void ComponentSystem::massConcentrationToMolarConcentration(const InputArrayType& mass_concentrations,
                                                            OutputArrayType& molar_concentrations) const {
    // Check dimensions using asserts (active only in Debug mode)
    assert(mass_concentrations.cols() == n_components &&
           "Number of columns in mass_concentrations must match number of components");
    assert(molar_concentrations.rows() == mass_concentrations.rows() &&
           molar_concentrations.cols() == mass_concentrations.cols() &&
           "Output array dimensions must match input array dimensions");

    // If molar mass is NaN, the result will also be NaN (IEEE 754 behavior)
    molar_concentrations.array() = mass_concentrations.array().rowwise() / molar_masses.array();
}

template <EigenArrayLike InputArrayType>
auto ComponentSystem::massConcentrationToMolarConcentration(const InputArrayType& mass_concentrations) const -> Array {
    Array result(mass_concentrations.rows(), mass_concentrations.cols());
    massConcentrationToMolarConcentration(mass_concentrations, result);
    return result;
}

// Template implementations for molarConcentrationToMassConcentration
template <EigenArrayLike InputArrayType, WritableEigenArrayLike OutputArrayType>
void ComponentSystem::molarConcentrationToMassConcentration(const InputArrayType& molar_concentrations,
                                                            OutputArrayType& mass_concentrations) const {
    // Check dimensions using asserts (active only in Debug mode)
    assert(molar_concentrations.cols() == n_components &&
           "Number of columns in molar_concentrations must match number of components");
    assert(mass_concentrations.rows() == molar_concentrations.rows() &&
           mass_concentrations.cols() == molar_concentrations.cols() &&
           "Output array dimensions must match input array dimensions");

    // If molar mass is NaN, the result will also be NaN (IEEE 754 behavior)
    mass_concentrations.array() = molar_concentrations.array().rowwise() * molar_masses.array();
}

template <EigenArrayLike InputArrayType>
auto ComponentSystem::molarConcentrationToMassConcentration(const InputArrayType& molar_concentrations) const -> Array {
    Array result(molar_concentrations.rows(), molar_concentrations.cols());
    molarConcentrationToMassConcentration(molar_concentrations, result);
    return result;
}

#endif  // COMPONENTSYSTEM_TPP
