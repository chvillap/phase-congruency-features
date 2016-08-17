/**
 * @file   main.cpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup phase_congruency
 * @ingroup    phase_congruency
 *
 * @copyright Copyright (c) 2016 Carlos Henrique Villa Pinto
 * @license GPL v2.0
 */

#include <cstdio>
#include <cstdlib>
#include <string>
#include "CmdLine.h"
#include "types.h"
#include "log_gabor_filter_bank.h"
#include "phase_congruency.h"


// Default argument values.
#define DEFAULT_NOISE_THRESHOLD "-1.0"
#define DEFAULT_NOISE_STD       "3.0" 
#define DEFAULT_SIGMOID_GAIN    "10.0"
#define DEFAULT_SIGMOID_CUTOFF  "0.5"


void show_usage(char *argv[])
{
    printf("----------------------------------------------------------------\n"
           "                        Phase Congruency                        \n"
           "----------------------------------------------------------------\n"
           " Basic usage: %s                                                \n"
           "              -i  input_image                                   \n"
           "              -fb filter_bank                                   \n"
           "              -o  output_filename_prefix                        \n"
           "----------------------------------------------------------------\n"
           " Options [description = default values]                         \n"
           "----------------------------------------------------------------\n"
           " Phase congruency parameters:                                   \n"
           "    -nt   [Noise energy suppression threshold = %s]             \n"
           "    -nstd [Noise standard deviation = %s]                       \n"
           "    -sig  [Sigmoid weighting parameters (gain cutoff) = %s %s]  \n"
           " Other:                                                         \n"
           "    -mask [Logical mask for the region of interest = NULL]      \n"
           "    -hann [Apply Hanning window function = true if given]       \n"
           "                                                                \n"
           " NOTE: for automatic noise threshold estimation, use -nt < 0    \n"
           "----------------------------------------------------------------\n\n",
        argv[0],
        DEFAULT_NOISE_THRESHOLD,
        DEFAULT_NOISE_STD,
        DEFAULT_SIGMOID_GAIN,
        DEFAULT_SIGMOID_CUTOFF
    );
}


int main(int argc, char *argv[])
{
    CCmdLine cmd;

    // Check for the help flag or for the lack of minimum required arguments.
    // Show the program usage if needed, and then exit.
    if (cmd.HasSwitch("-h") || cmd.SplitLine(argc, argv) < 1) {
        show_usage(argv);
        return EXIT_FAILURE;
    }

    // ------------------------------------------------------------------------
    // Getting the inputs:

    std::string filename_input;
    std::string filename_filter_bank;
    std::string filename_output_prefix;
    try
    {
        filename_input         = cmd.GetArgument("-i", 0);
        filename_filter_bank   = cmd.GetArgument("-fb", 0);
        filename_output_prefix = cmd.GetArgument("-o", 0);
    }
    catch (...)
    {
        show_usage(argv);
        return EXIT_FAILURE;
    }

    float noise_threshold = static_cast<float>(atof(
        cmd.GetSafeArgument("-nt", 0, DEFAULT_NOISE_THRESHOLD).c_str()));
    float noise_std = static_cast<float>(atof(
        cmd.GetSafeArgument("-nstd", 0, DEFAULT_NOISE_STD).c_str()));
    float sigmoid_gain = static_cast<float>(atof(
        cmd.GetSafeArgument("-sig", 0, DEFAULT_SIGMOID_GAIN).c_str()));
    float sigmoid_cutoff = static_cast<float>(atof(
        cmd.GetSafeArgument("-sig", 1, DEFAULT_SIGMOID_CUTOFF).c_str()));

    std::string filename_mask = cmd.GetSafeArgument("-mask", 0, "");
    bool hanning_flag = cmd.HasSwitch("-hann");

    /*
     * TODO:
     * Use some library to read 2D/3D images here.
     * And then get the image size.
     */
    float *input_image = NULL;
    bip::triple<size_t> sizes;

    bool *input_mask = NULL;
    /*
     * TODO:
     * Use some library to read 2D/3D images here.
     * And then compare the mask size to the input image size.
     */

    bip::log_gabor_filter_bank *lgfb =
        bip::log_gabor_filter_bank::read_parameters(filename_filter_bank);

    // ------------------------------------------------------------------------
    // Computing the phase congruency and saving the results:

    if (hanning_flag) {
        for (size_t i = 0, z = 0; z < sizes[2]; ++z)
            for (size_t y = 0; y < sizes[1]; ++y)
                for (size_t x = 0; x < sizes[0]; ++x, ++i) {
                    bip::triple<size_t> xyz = {x, y, z};
                    input_image[i] *= bip::hanning(xyz, sizes);
                }
    }

    bip::phase_congruency pc(filename_output_prefix, input_image, lgfb, sizes,
                             input_mask, noise_threshold, noise_std,
                             sigmoid_gain, sigmoid_cutoff);
    pc.compute();

    // ------------------------------------------------------------------------
    // Cleaning up:

    if (input_mask != NULL)
        delete[] input_mask;
    delete[] input_image;
    delete lgfb;

    return EXIT_SUCCESS;
}
