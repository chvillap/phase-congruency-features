/**
 * @file   math_functions.h
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup utils
 * @ingroup    utils
 *
 * @copyright Copyright (c) 2016 Carlos Henrique Villa Pinto
 * @license GPL v2.0
 */

#ifndef MATH_FUNCTIONS_H
#define MATH_FUNCTIONS_H

#define _USE_MATH_DEFINES
#include <cstddef>
#include <cfloat>
#include <cmath>
#include <algorithm>
#include <vector>
#include <fftw3.h>
#include "assert2.h"
#include "types.h"


namespace bip
{


/**
 * @var EPSILON
 *
 * @brief A very small number commonly used to avoid division by zero.
 */
extern const float EPSILON;

/**
 * @fn fp_equal
 *
 * @brief Equality test for floating-point numbers.
 * 
 * @param[in] n1       First operand.
 * @param[in] n2       Second operand.
 * @param[in] maxerror Maximum error below which n1 and n2 are considered equal.
 *
 * @returns True if n1 and n2 are equal or close enough; false otherwise.
 *
 * @attention This function should always be used in place of the == operator
 * for floating-point numbers.
 */
bool fp_equal(float n1, float n2, float maxerror = EPSILON);

/**
 * @fn sqr
 *
 * @brief Computes the square of a real number.
 * 
 * @param[in] n Input number.
 *
 * @returns The square n^2 of the input number.
 */
inline
float sqr(float n) {
    return n * n;
}

/**
 * @fn sign
 *
 * @brief Sign function of a real number.
 * 
 * @param[in] n Input number.
 *
 * @returns -1 for negative numbers, 1 for positive numbers, 0 for zero.
 */
inline
int sign(float n) {
    return (n > 0) - (n < 0);
}

/**
 * @fn posmod
 *
 * @brief Positive modulus (remainder of two real numbers).
 * 
 * @param[in] n First operand (dividend).
 * @param[in] m Second operand (divisor).
 *
 * @returns The positive remainder r = n - m * floor(n/m).
 */
inline
float posmod(float n, float m) {
    const float q = n / m;
    return (q - floor(q)) * m;
}

/**
 * @fn rad2deg
 *
 * @brief Converts radians to degrees.
 * 
 * @param[in] rad A measure in radians.
 *
 * @returns The same measure in degrees.
 */
inline
float rad2deg(float rad) {
    return rad / M_PI * 180;
}

/**
 * @fn deg2rad
 *
 * @brief Converts degrees to radians.
 *
 * @param[in] deg A measure in degrees.
 *
 * @returns The same measure in radians.
 */
inline
float deg2rad(float deg) {
    return deg / 180 * M_PI;
}

/**
 * @fn cart2sph
 *
 * @brief Converts cartesian coordinates to spherical coordinates.
 * 
 * @param[in] cart A position in Cartesian coordinates (x, y, z).
 *
 * @returns The same position in spherical coordinates (rho, phi, theta).
 *
 * @warning The spherical elevation angle (theta) goes from -pi to pi.
 */
triple<float> cart2sph(triple<float> cart);

/**
 * @fn sph2cart
 *
 * @brief Converts spherical coordinates to cartesian coordinates.
 * 
 * @param[in] sph Spherical coordinates (rho, phi, theta).
 *
 * @returns The same position in cartesian coordinates (x, y, z).
 *
 * @warning The spherical elevation angle goes from -pi to pi.
 */
triple<float> sph2cart(triple<float> sph);

/**
 * @fn FFT
 *
 * @brief Computes the forward/backward Fast Fourier Transform (FFT) on a 2D/3D
 *        matrix.
 * 
 * @param[in] M        Input matrix.
 * @param[in] sizes    Number of matrix elements by dimension.
 * @param[in] backward Tells whether it should compute the backward FFT.
 *
 * @warning The FFT algorithm works faster for matrices whose dimensions are
 * powers of two.
 *
 * @see http://www.fftw.org/
 */
void FFT(fftwf_complex *M, triple<size_t> sizes, bool backward = false);

/**
 * @fn hanning
 *
 * @brief Computes the Hanning window function value at a given position.
 * 
 * @param[in] xyz   Position in cartesian coordinates.
 * @param[in] sizes Domain size for each dimension.
 *
 * @returns A real value in the [0-1] range.
 */
float hanning(triple<size_t> xyz, triple<size_t> sizes);

/**
 * @fn median
 *
 * @brief Computes the median value in an array.
 * 
 * @param[in] array Input array.
 * @param[in] size  Number of elements in the array.
 *
 * @returns The median value.
 */
float median(float *array, size_t size);

/**
 * @fn normalize_min_max
 *
 * @brief Normalizes array values to the [new_min, new_max] range.
 * 
 * @param[in] array   Input array.
 * @param[in] size    Number of elements in the array.
 * @param[in] new_min Target minimum value after normalization.
 * @param[in] new_max Target maximum value after normalization.
 */
void normalize_min_max(float *array,
                       size_t size,
                       float  new_min,
                       float  new_max);


}

#endif
