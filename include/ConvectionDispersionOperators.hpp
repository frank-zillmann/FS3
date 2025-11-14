#ifndef CONVECTION_DISPERSION_OPERATORS_HPP
#define CONVECTION_DISPERSION_OPERATORS_HPP

#include "EigenDataTypes.hpp"

// CONVECTION:
// The following header-only templated kernels apply 1D finite-volume
// convection (upwind / WENO3) and dispersion (central differences) on
// uniform or non-uniform cells. They operate elementwise on Eigen arrays
// representing masses per cell and accumulate into dm_dt.

// Templated versions that work with any Eigen array-like types
template <EigenArrayLike MArrayType, WritableEigenArrayLike DmDtArrayType>
void internalUniformUpwind(const MArrayType& m, DmDtArrayType& dm_dt, realtype flowRate, realtype cell_volume);

template <EigenArrayLike MArrayType, WritableEigenArrayLike DmDtArrayType, EigenArrayLike V_l_ArrayType>
void internalNonUniformUpwind(const MArrayType& m, DmDtArrayType& dm_dt, realtype flowRate, const V_l_ArrayType& V_l);

template <EigenArrayLike MArrayType, WritableEigenArrayLike DmDtArrayType>
void internalUniformWENO3(const MArrayType& m, DmDtArrayType& dm_dt, realtype flowRate, realtype cell_volume);

template <EigenArrayLike MFromArrayType, WritableEigenArrayLike DmDtFromArrayType, EigenArrayLike MToArrayType, WritableEigenArrayLike DmDtToArrayType>
void connectionUpwind(const MFromArrayType& m_from,
                      DmDtFromArrayType& dm_dt_from,
                      realtype cell_volume_from,
                      const MToArrayType& m_to,
                      DmDtToArrayType& dm_dt_to,
                      realtype flowRate);

// DISPERSION:
template <EigenArrayLike MArrayType, WritableEigenArrayLike DmDtArrayType>
void internalUniformDispersion(const MArrayType& m,
                               DmDtArrayType& dm_dt,
                               realtype dispersion_coefficient,
                               realtype cell_distance);

template <EigenArrayLike MArrayType, WritableEigenArrayLike DmDtArrayType, EigenArrayLike VolumeArrayType>
void internalNonUniformDispersion(const MArrayType& m,
                                  DmDtArrayType& dm_dt,
                                  const VolumeArrayType& V_l,
                                  realtype dispersion_coefficient,
                                  realtype cell_distance);

// Include template implementations
#include "ConvectionDispersionOperators.tpp"

#endif  // CONVECTION_DISPERSION_OPERATORS_HPP