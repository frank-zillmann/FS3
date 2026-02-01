/**
 * @file NumpyIO.hpp
 * @brief Wrapper functions for saving data to NumPy .npy and .npz files
 *
 * These functions provide a simple interface to cnpy's functionality,
 * allowing FS3 users to save simulation data in NumPy-compatible formats.
 */

#ifndef FS3_NUMPY_IO_HPP
#define FS3_NUMPY_IO_HPP

#include <string>
#include <vector>

namespace FS3 {

/**
 * @brief Save data to a .npy file (NumPy array format)
 *
 * @tparam T Data type (double, float, int, etc.)
 * @param filename Path to the output file (should end with .npy)
 * @param data Pointer to the data array
 * @param shape Shape of the array (e.g., {rows, cols} for 2D)
 * @param mode "w" to overwrite, "a" to append along first axis
 *
 * @example
 * std::vector<double> data = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
 * FS3::npy_save("output.npy", data.data(), {2, 3}, "w");
 */
template <typename T>
void npy_save(const std::string& filename,
              const T* data,
              const std::vector<size_t>& shape,
              const std::string& mode = "w");

/**
 * @brief Save data to a .npz file (compressed NumPy archive)
 *
 * @tparam T Data type (double, float, int, etc.)
 * @param zipname Path to the output file (should end with .npz)
 * @param varname Name of the variable inside the archive
 * @param data Pointer to the data array
 * @param shape Shape of the array (e.g., {rows, cols} for 2D)
 * @param mode "w" to overwrite/create, "a" to append variable to archive
 *
 * @example
 * std::vector<double> timestamps = {0.0, 1.0, 2.0};
 * std::vector<double> values = {10.0, 20.0, 30.0};
 * FS3::npz_save("results.npz", "timestamps", timestamps.data(), {3}, "w");
 * FS3::npz_save("results.npz", "values", values.data(), {3}, "a");
 */
template <typename T>
void npz_save(const std::string& zipname,
              const std::string& varname,
              const T* data,
              const std::vector<size_t>& shape,
              const std::string& mode = "w");

}  // namespace FS3

#endif  // FS3_NUMPY_IO_HPP
