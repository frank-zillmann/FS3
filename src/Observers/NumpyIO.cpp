/**
 * @file NumpyIO.cpp
 * @brief Implementation of NumPy I/O wrapper functions
 */

#include "Observers/NumpyIO.hpp"

#include "cnpy.h"

namespace FS3 {

// Explicit instantiations for common types
template <typename T>
void npy_save(const std::string& filename, const T* data, const std::vector<size_t>& shape, const std::string& mode) {
    cnpy::npy_save(filename, data, shape, mode);
}

template <typename T>
void npz_save(const std::string& zipname,
              const std::string& varname,
              const T* data,
              const std::vector<size_t>& shape,
              const std::string& mode) {
    cnpy::npz_save(zipname, varname, data, shape, mode);
}

// Explicit template instantiations for common numeric types
template void npy_save<double>(const std::string&, const double*, const std::vector<size_t>&, const std::string&);
template void npy_save<float>(const std::string&, const float*, const std::vector<size_t>&, const std::string&);
template void npy_save<int>(const std::string&, const int*, const std::vector<size_t>&, const std::string&);
template void npy_save<long>(const std::string&, const long*, const std::vector<size_t>&, const std::string&);
template void npy_save<size_t>(const std::string&, const size_t*, const std::vector<size_t>&, const std::string&);

template void npz_save<double>(const std::string&,
                               const std::string&,
                               const double*,
                               const std::vector<size_t>&,
                               const std::string&);
template void npz_save<float>(const std::string&,
                              const std::string&,
                              const float*,
                              const std::vector<size_t>&,
                              const std::string&);
template void npz_save<int>(const std::string&, const std::string&, const int*, const std::vector<size_t>&, const std::string&);
template void npz_save<long>(const std::string&, const std::string&, const long*, const std::vector<size_t>&, const std::string&);
template void npz_save<size_t>(const std::string&,
                               const std::string&,
                               const size_t*,
                               const std::vector<size_t>&,
                               const std::string&);

}  // namespace FS3
