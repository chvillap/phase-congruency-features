/**
 * @file   assert2.cpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup utils
 * @ingroup    utils
 *
 * @copyright Copyright (c) 2016 Carlos Henrique Villa Pinto
 * @license GPL v2.0
 */

#include "assert2.h"


namespace bip
{
namespace debug
{


void
_assert2(bool passed, const char *assertion, const char *file, long line)
{
    if (!passed) {
        std::stringstream ss;
        ss << "Failed assertion " << assertion
           << " in file " << file
           << " at line " << line;

        std::string message(ss.str());
        std::cerr << message << std::endl;

        throw std::logic_error(message);
    }
}


}
}
