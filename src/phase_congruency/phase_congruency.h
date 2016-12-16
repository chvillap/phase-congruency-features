/**
 * @file   phase_congruency.h
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup phase_congruency
 * @ingroup    phase_congruency
 *
 * @copyright Copyright (c) 2016 Carlos Henrique Villa Pinto
 * @license GPL v2.0
 */

#ifndef PHASE_CONGRUENCY_H
#define PHASE_CONGRUENCY_H

#define DEFAULT_FFTW_NUMBER_OF_THREADS 16
#define _USE_MATH_DEFINES
#include <cstddef>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>
#include <new>
#include <fftw3.h>
#include "assert2.h"
#include "types.h"
#include "math_functions.h"
#include "image_io.h"
#include "log_gabor_filter_bank.h"

#define PHASE_CONGRUENCY_VERBOSE_ON
// #define PHASE_CONGRUENCY_DEBUG_ON

namespace bip
{


/**
 * @class phase_congruency phase_congruency.h
 *
 * @brief The phase congruency model for invariant 2D/3D image feature
 * detection.
 *
 * @see Kovesi, P., 2000. Phase congruency:
 * A low-level image invariant. Psychological Research 64, 136-148.
 */
class phase_congruency
{
public:
    /**
     * @brief Constructor method.
     *
     * @param[in] filename_prefix Prefix of the output file names.
     * @param[in] input_image     Input image data.
     * @param[in] filter_bank     Bank of log-Gabor filters.
     * @param[in] sizes           Image and filter size in each dimension.
     * @param[in] input_mask      Logical mask of the region of interest.
     * @param[in] noise_threshold Noise energy suppression threshold.
     * @param[in] noise_std       Standard deviation of the noise model.
     * @param[in] sigmoid_gain    Sigmoidal response weighting gain.
     * @param[in] sigmoid_cutoff  Sigmoidal response weighting cut-off.
     */
    phase_congruency(std::string            filename_prefix,
                     float                 *input_image,
                     log_gabor_filter_bank *filter_bank,
                     triple<size_t>         sizes,
                     bool                  *input_mask      = NULL,
                     float                  noise_threshold = -1.0,
                     float                  noise_std       = 3.0,
                     float                  sigmoid_gain    = 10.0,
                     float                  sigmoid_cutoff  = 0.5);

    /**
     * @brief Destructor method.
     */
    ~phase_congruency();

    /**
     * @brief Prints a string representation of the object.
     */
    std::ostream& print(std::ostream &os) const;

    /**
     * @brief Computes the phase congruency technique on the input image.
     */
    void compute();

    /** @brief Getter for m_filename_prefix. */
    std::string get_filename_prefix() const {
        return m_filename_prefix;
    }

    /** @brief Getter for m_filter_bank. */
    log_gabor_filter_bank* get_filter_bank() const {
        return m_filter_bank;
    }

    /** @brief Getter for m_input_image. */
    float* get_input_image() const {
        return m_input_image;
    }

    /** @brief Getter for m_input_mask. */
    bool* get_input_mask() const {
        return m_input_mask;
    }

    /** @brief Getter for m_sizes. */
    triple<size_t> get_sizes() const {
        return m_sizes;
    }

    /** @brief Getter for m_noise_threshold. */
    float get_noise_threshold() const {
        return m_noise_threshold;
    }

    /** @brief Getter for m_noise_std. */
    float get_noise_std() const {
        return m_noise_std;
    }

    /** @brief Getter for m_sigmoid_gain. */
    float get_sigmoid_gain() const {
        return m_sigmoid_gain;
    }

    /** @brief Getter for m_sigmoid_cutoff. */
    float get_sigmoid_cutoff() const {
        return m_sigmoid_cutoff;
    }

private:
    /**
     * @brief Computes the DC-shifted forward FFT of the input image.
     *
     * @param[out] target Preallocated complex image to compute.
     *
     * @return The resulting Fourier transform as a complex image.
     */
    void compute_shifted_FFT(fftwf_complex *f_target);

    /**
     * @brief Applies a single log-Gabor filter to a complex image in the
     * frequency domain.
     *
     * @param[out] f_output  Preallocated output complex image.
     * @param[in]  f_input   Preallocated input complex image.
     * @param[in]  scale     Filter scale parameter.
     * @param[in]  azimuth   Filter azimuth parameter.
     * @param[in]  elevation Filter elevation parameter.
     */
    void compute_filtering(fftwf_complex *f_output,
                           fftwf_complex *f_input,
                           size_t         scale,
                           size_t         azimuth,
                           size_t         elevation);

    /**
     * @brief Automatically estimates the noise energy threshold from the
     * amplitudes of filter responses.
     *
     * @param[in] sum_amplitude Array of sums of filter response amplitudes
     * accumulated over all the scales of the bank of filters.
     *
     * @return The estimated noise energy threshold.
     */
    float estimate_noise_threshold(float *sum_amplitude);

    /**
     * @brief Applies an energy weighting function to suppress responses that
     * locally occur only to a few frequency components.
     *
     * @param[in] energy        Local energy value.
     * @param[in] sum_amplitude Local summed filter response amplitude.
     * @param[in] max_amplitude Local maximum filter response amplitude.
     *
     * @return The weighted local energy value, as a positive real number.
     */
    float apply_energy_weighting(float energy,
                                 float sum_amplitude,
                                 float max_amplitude);

    /** @brief Setter for m_filename_prefix. */
    void set_filename_prefix(std::string filename_prefix) {
        debug::assert2(!filename_prefix.empty());
        m_filename_prefix = filename_prefix;
    }

    /** @brief Setter for m_filter_bank. */
    void set_filter_bank(log_gabor_filter_bank *filter_bank) {
        debug::assert2(filter_bank != NULL);
        m_filter_bank = filter_bank;
    }

    /** @brief Setter for m_input_image. */
    void set_input_image(float *input_image) {
        debug::assert2(input_image != NULL);
        m_input_image = input_image;
    }

    /** @brief Setter for m_input_mask. */
    void set_input_mask(bool *input_mask) {
        m_input_mask = input_mask;
    }

    /** @brief Setter for m_sizes. */
    void set_sizes(triple<size_t> sizes) {
        m_sizes = sizes;
    }

    /** @brief Setter for m_noise_threshold. */
    void set_noise_threshold(float noise_threshold) {
        m_noise_threshold = noise_threshold;
    }

    /** @brief Setter for m_noise_std. */
    void set_noise_std(float noise_std) {
        debug::assert2(noise_std > 0.0);
        m_noise_std = noise_std;
    }

    /** @brief Setter for m_sigmoid_gain. */
    void set_sigmoid_gain(float sigmoid_gain) {
        m_sigmoid_gain = sigmoid_gain;
    }

    /** @brief Setter for m_sigmoid_cutoff. */
    void set_sigmoid_cutoff(float sigmoid_cutoff) {
        debug::assert2(sigmoid_cutoff >= 0.0 &&
                       sigmoid_cutoff <= 1.0);
        m_sigmoid_cutoff = sigmoid_cutoff;
    }

    /** @attention Purposely not implemented. */
    phase_congruency(const phase_congruency &other);

    /** @attention Purposely not implemented. */
    phase_congruency& operator=(const phase_congruency &rhs);

private:
    std::string            m_filename_prefix;
    log_gabor_filter_bank *m_filter_bank;
    triple<size_t>         m_sizes;
    float                 *m_input_image;
    bool                  *m_input_mask;
    float                  m_noise_threshold;
    float                  m_noise_std;
    float                  m_sigmoid_gain;
    float                  m_sigmoid_cutoff;
};


// ------------------------------------------------------------------------------------------------

/**
 * @brief Output stream operator for bip::phase_congruency.
 */
std::ostream& operator<<(std::ostream &os, const phase_congruency &pc);


}

#endif
