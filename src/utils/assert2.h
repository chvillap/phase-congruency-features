/**
 * @file   assert.h
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup utils
 * @ingroup    utils
 *
 * @copyright Copyright (c) 2016 Carlos Henrique Villa Pinto
 * @license GPL v2.0
 */

#ifndef ASSERT2_H
#define ASSERT2_H

#include <cstddef>
#include <cassert>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <string>


namespace bip
{
namespace debug
{


/**
 * @fn _assert2
 *
 * @brief Internal function to display an assertion failed message.
 *
 * @param[in] passed    Tells whether the assertion has passed or not.
 * @param[in] assertion Evaluated assertion.
 * @param[in] file      Source file name.
 * @param[in] line      Line of code in the source file.
 *
 * @attention This function is automatically called by the bip::debug::assert
 *            macro and therefore it should not be directly called by the user.
 *
 * @see assert2
 */
void _assert2(bool passed, const char *assertion, const char *file, long line);

/**
 * @def assert2
 *
 * @brief Macro for assertion checking that works for any build type.
 *
 * @param[in] expr Assertion expression.
 *
 * @see _assert2
 */
// #ifdef NDEBUG
    // #define assert2(expr) _assert2(true, "", "", 0)
// #else
    #define assert2(expr) _assert2(expr, #expr, __FILE__, __LINE__)
// #endif


}
}

#endif
