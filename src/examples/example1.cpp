/**
 * @file   example1.cpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup examples
 * @ingroup    examples
 *
 * @copyright Copyright (c) 2016 Carlos Henrique Villa Pinto
 * @license GPL v2.0
 */

#include <iostream>
#include <string>
#include "log_gabor_filter_bank.h"
#include "image_io.h"
#include "types.h"


int main(int argc, char *argv[])
{
    // Generate a bank of 2D log-Gabor filters.
    try {
        bip::triple<size_t> size = {128, 128, 1};
        bip::log_gabor_filter_bank lgbf_1(
            "log_gabor", // Filename prefix.
            size,           // Filter size (z=1 for 2D).
            3,              // Scales.
            6,              // Azimuths.
            1,              // Elevations (1 for 2D).
            1./3,           // Max central frequency.
            2.1,            // Multiplicative factor.
            0.55,           // Frequency spread ratio.
            1.2,            // Angular spread ratio.
            15,             // Butterworth order.
            0.45,           // Butterworth cutoff.
            false           // Uniform sampling?
        );
        bip::log_gabor_filter_bank::write_parameters(lgbf_1);
        lgbf_1.compute();

        std::cout << lgbf_1 << std::endl;

        std::string filename(
            "log_gabor_#_128_128_001_03_06_01_033_209_055_120_0.bof");

        bip::log_gabor_filter_bank *lgbf_2 =
            bip::log_gabor_filter_bank::read_parameters(filename);

        std::cout << *lgbf_2 << std::endl;

        delete lgbf_2;

    } catch (const char *exception) {
        std::cerr << exception << std::endl;
    }

    return EXIT_SUCCESS;
}
