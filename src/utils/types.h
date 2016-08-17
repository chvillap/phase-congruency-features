/**
 * @file   types.h
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup utils
 * @ingroup    utils
 *
 * @copyright Copyright (c) 2016 Carlos Henrique Villa Pinto
 * @license GPL v2.0
 */

#ifndef TYPES_H
#define TYPES_H

#include <array>
#include <iostream>


namespace bip
{


/**
 * @def triple
 *
 * @brief A generic array of size 3.
 * 
 * @warning This is just an alias for the std::array class, so C++11 must be
 *          supported by the compiler.
 */
template <class T>
using triple = std::array<T, 3>;

// ----------------------------------------------------------------------------

/**
 * @brief Output stream operator for bip::triple.
 */
template <class T>
std::ostream&
operator<<(std::ostream &os, const triple<T> &t)
{
    return os << "(" << t[0] << ", " << t[1] << ", " << t[2] << ")";
}


}

#endif
