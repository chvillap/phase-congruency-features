/**
 * @file   example2.cpp
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
#include "types.h"
#include "image_io.h"
#include "log_gabor_filter_bank.h"
#include "phase_congruency.h"


int main(int argc, char *argv[])
{
    typedef itk::Image<float, 3> TImage;

    // Compute the phase congruency for a 2D image.
    // First it creates a bank of 2D log-Gabor filters (you can skip this
    // if the filters already exist in disk).
    try {
        bip::triple<size_t> size = {256, 256, 1};

        bip::log_gabor_filter_bank lgbf(
            "log_gabor", // Filename prefix.
            size,        // Filter size (z=1 for 2D).
            3,           // Scales.
            6,           // Azimuths.
            1,           // Elevations (1 for 2D).
            1./3,        // Max central frequency.
            2.1,         // Multiplicative factor.
            0.55,        // Frequency spread ratio.
            1.2,         // Angular spread ratio.
            15,          // Butterworth order.
            0.45,        // Butterworth cutoff.
            false        // Uniform sampling?
        );
        bip::log_gabor_filter_bank::write_parameters(lgbf);
        lgbf.compute();

        TImage::Pointer itk_image = bip::io::read_image<float, 3>(
            "data/cameraman.tif");
        float *image = bip::io::image2array<float, 3>(itk_image);

        bip::phase_congruency pc(
            "cameraman", // Filename prefix.
            image,       // Input image.
            &lgbf,       // Bank of log-gabor filters.
            size,        // Image size (z=1 for 2D).
            NULL,        // Input mask (NULL for no mask).
            -1.0,        // Noise energy threshold (< 0 for auto estimation).
            2.0,         // Noise standard deviation.
            10,          // Sigmoid weighting gain.
            0.5          // Sigmoid weighting cutoff.
        );
        pc.compute();

        delete[] image;

    } catch (const char *exception) {
        std::cerr << exception << std::endl;
    }

    return EXIT_SUCCESS;
}
