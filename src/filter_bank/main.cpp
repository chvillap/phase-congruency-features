/**
 * @file   main.cpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup filter_bank
 * @ingroup    filter_bank
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


// Default argument values.
const char *DEFAULT_IMAGE_SIZE_X    = "128";
const char *DEFAULT_IMAGE_SIZE_Y    = "128";
const char *DEFAULT_IMAGE_SIZE_Z    = "128";
const char *DEFAULT_NUM_SCALES      = "4";
const char *DEFAULT_NUM_AZIMUTHS    = "6";
const char *DEFAULT_NUM_ELEVATIONS  = "3";
const char *DEFAULT_MAX_FREQUENCY   = "0.333333";
const char *DEFAULT_MULT_FACTOR     = "2.1";
const char *DEFAULT_FREQUENCY_RATIO = "0.55";
const char *DEFAULT_ANGULAR_RATIO   = "1.2";
const char *DEFAULT_LOWPASS_ORDER   = "15";
const char *DEFAULT_LOWPASS_CUTOFF  = "0.45";


void show_usage(char *argv[])
{
    printf("----------------------------------------------------------------\n"
           "                      Log-Gabor Filter Bank                     \n"
           "----------------------------------------------------------------\n"
           " Basic usage: %s                                                \n"
           "              -o output-filename-prefix                         \n"
           "----------------------------------------------------------------\n"
           " Options [description = default values]                         \n"
           "----------------------------------------------------------------\n"
           " Image parameters:                                              \n"
           "    -sz   [Image sizes (cols rows slices) = %s %s %s]           \n"
           " Filter bank parameters:                                        \n"
           "    -ns   [Number of scales = %s]                               \n"
           "    -na   [Number of azimuth angles = %s]                       \n"
           "    -ne   [Number of elevation angles = %s]                     \n"
           "    -maxf [Maximum central frequency = %s]                      \n"
           "    -mulf [Multiplicative factor = %s]                          \n"
           "    -fr   [Frequency bandwidth ratio = %s]                      \n"
           "    -ar   [Angular spread ratio = %s]                           \n"
           " Butterworth lowpass filter parameters:                         \n"
           "    -lpo  [Lowpass filter order = %s]                           \n"
           "    -lpc  [Lowpass filter cutoff = %s]                          \n"
           " Other:                                                         \n"
           "    -us   [Uniform sampling = true if given]                    \n"
           "                                                                \n"
           " NOTE: for 2D images, use -sz cols rows 1 and -ne 1             \n"
           "----------------------------------------------------------------\n\n",
        argv[0],
        DEFAULT_IMAGE_SIZE_X,
        DEFAULT_IMAGE_SIZE_Y,
        DEFAULT_IMAGE_SIZE_Z,
        DEFAULT_NUM_SCALES,
        DEFAULT_NUM_AZIMUTHS,
        DEFAULT_NUM_ELEVATIONS,
        DEFAULT_MAX_FREQUENCY,
        DEFAULT_MULT_FACTOR,
        DEFAULT_FREQUENCY_RATIO,
        DEFAULT_ANGULAR_RATIO,
        DEFAULT_LOWPASS_ORDER,
        DEFAULT_LOWPASS_CUTOFF
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

    std::string filename_prefix;
    try
    {
        filename_prefix = cmd.GetArgument("-o", 0);
    }
    catch (...)
    {
        show_usage(argv);
        return EXIT_FAILURE;
    }

    size_t size_x = static_cast<size_t>(atol(
        cmd.GetSafeArgument("-sz", 0, DEFAULT_IMAGE_SIZE_X).c_str()));
    size_t size_y = static_cast<size_t>(atol(
        cmd.GetSafeArgument("-sz", 1, DEFAULT_IMAGE_SIZE_Y).c_str()));
    size_t size_z = static_cast<size_t>(atol(
        cmd.GetSafeArgument("-sz", 2, DEFAULT_IMAGE_SIZE_Z).c_str()));

    size_t num_scales = static_cast<size_t>(atol(
        cmd.GetSafeArgument("-ns", 0, DEFAULT_NUM_SCALES).c_str()));
    size_t num_azimuths = static_cast<size_t>(atol(
        cmd.GetSafeArgument("-na", 0, DEFAULT_NUM_AZIMUTHS).c_str()));
    size_t num_elevations = static_cast<size_t>(atol(
        cmd.GetSafeArgument("-ne", 0, DEFAULT_NUM_ELEVATIONS).c_str()));
    
    float max_frequency = static_cast<float>(atof(
        cmd.GetSafeArgument("-maxf", 0, DEFAULT_MAX_FREQUENCY).c_str()));
    float mult_factor = static_cast<float>(atof(
        cmd.GetSafeArgument("-mulf", 0, DEFAULT_MULT_FACTOR).c_str()));

    float frequency_ratio = static_cast<float>(atof(
        cmd.GetSafeArgument("-fr", 0, DEFAULT_FREQUENCY_RATIO).c_str()));
    float angular_ratio = static_cast<float>(atof(
        cmd.GetSafeArgument("-ar", 0, DEFAULT_ANGULAR_RATIO).c_str()));
    
    float lowpass_order = static_cast<float>(atof(
        cmd.GetSafeArgument("-lpo", 0, DEFAULT_LOWPASS_ORDER).c_str()));
    float lowpass_cutoff = static_cast<float>(atof(
        cmd.GetSafeArgument("-lpc", 0, DEFAULT_LOWPASS_CUTOFF).c_str()));
    
    bool uniform_sampling = cmd.HasSwitch("-us");

    // ------------------------------------------------------------------------
    // Computing/saving the log-Gabor filter bank:

    bip::triple<size_t> sizes = {size_x, size_y, size_z};
    
    bip::log_gabor_filter_bank lgfb(filename_prefix, sizes, num_scales,
                                    num_azimuths, num_elevations,
                                    max_frequency, mult_factor,
                                    frequency_ratio, angular_ratio,
                                    lowpass_order, lowpass_cutoff,
                                    uniform_sampling);
    
    bip::log_gabor_filter_bank::write_parameters(lgfb);
    lgfb.compute();

    return EXIT_SUCCESS;
}
