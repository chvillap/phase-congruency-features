/**
 * @file   log_gabor_filter_bank.h
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup filter_bank
 * @ingroup    filter_bank
 *
 * @copyright Copyright (c) 2016 Carlos Henrique Villa Pinto
 * @license GPL v2.0
 */

#ifndef LOG_GABOR_FILTER_BANK_H
#define LOG_GABOR_FILTER_BANK_H

#include <cstddef>
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "assert2.h"
#include "types.h"
#include "math_functions.h"

#define LOG_GABOR_FILTER_BANK_VERBOSE_ON
#define LOG_GABOR_FILTER_BANK_DEBUG_ON


namespace bip
{


/**
 * @class log_gabor_filter_bank log_gabor_filter_bank.h
 *
 * @brief A bank of 2D/3D log-Gabor filters generated in the frequency domain,
 * which can be easily saved as files for further reuse.
 * 
 * @attention Your filename prefix must not include the "_#_" sequence.
 *
 * @see https://en.wikipedia.org/wiki/Log_Gabor_filter
 * @see http://www.peterkovesi.com/matlabfns/PhaseCongruency/Docs/convexpl.html
 */
class log_gabor_filter_bank
{
public:
    /**
     * @brief Constructor method.
     * 
     * @param[in] filename_prefix  Prefix of the filter file names.
     * @param[in] sizes            Filter size in each dimension.
     * @param[in] num_scales       Number of filter scales.
     * @param[in] num_azimuths     Number of filter azimuth angles.
     * @param[in] num_elevations   Number of filter elevation angles.
     * @param[in] max_frequency    Maximum central frequency.
     * @param[in] mult_factor      Multiplicative factor of filter frequencies.
     * @param[in] frequency_ratio  Frequency spread ratio.
     * @param[in] angular_ratio    Angular spread ratio.
     * @param[in] lowpass_order    Butterworth lowpass filter order.
     * @param[in] lowpass_cutoff   Butterworth lowpass filter cut-off.
     * @param[in] uniform_sampling Defines the filter sampling approach.
     */
    log_gabor_filter_bank(std::string    filename_prefix,
                          triple<size_t> sizes,
                          size_t         num_scales       = 4,
                          size_t         num_azimuths     = 6,
                          size_t         num_elevations   = 3,
                          float          max_frequency    = 1./3,
                          float          mult_factor      = 2.1,
                          float          frequency_ratio  = 0.55,
                          float          angular_ratio    = 1.2,
                          float          lowpass_order    = 15.0,
                          float          lowpass_cutoff   = 0.45,
                          bool           uniform_sampling = false);

    /**
     * @brief Destructor method.
     */
    ~log_gabor_filter_bank();

    /**
     * @brief Prints a string representation of the object.
     */
    std::ostream& print(std::ostream &os) const;

    /**
     * @brief Get a specific filter of the bank from file.
     * 
     * @param[in] scale     Filter scale parameter.
     * @param[in] azimuth   Filter azimuth parameter.
     * @param[in] elevation Filter elevation parameter.
     *
     * @returns The filter data array.
     */
    float* get_filter(size_t scale,
                      size_t azimuth = 0,
                      size_t elevation = 0);

    /**
     * @brief Computes (generates and saves) all the filters of the bank.
     */
    void compute();

    /**
     * @brief Reads the parameters of the filter bank from a header file.
     * 
     * @param[in] filename Header file name.
     *
     * @returns An instance of the log-Gabor filter bank object.
     *
     * @attention The filename must not include the "_#_" sequence.
     */
    static log_gabor_filter_bank* read_parameters(std::string filename);

    /**
     * @brief Writes the parameters of the filter bank into a header file.
     * 
     * @param[in] bof Reference to the log-Gabor filter bank object.
     */
    static void write_parameters(log_gabor_filter_bank &bof);

    /**
     * @brief Calculates the total number of filter orientations according to
     * the numbers of azimuth and elevation angles and the filter sampling.
     * 
     * @returns The total number of filter orientations.
     */
    size_t get_num_orientations() const {
        size_t result = 0;
        for (size_t e = 0; e < m_num_elevations; ++e)
            result += m_num_azimuths_per_elevation[e];
        return result;
    }

    /** @brief Getter for m_filename_prefix. */
    std::string get_filename_prefix() const {
        return m_filename_prefix;
    }

    /** @brief Getter for m_sizes. */
    triple<size_t> get_sizes() const {
        return m_sizes;
    }

    /** @brief Getter for m_num_scales. */
    size_t get_num_scales() const {
        return m_num_scales;
    }

    /** @brief Getter for m_num_azimuths. */
    size_t get_num_azimuths() const {
        return m_num_azimuths;
    }

    /** @brief Getter for m_num_azimuths_per_elevation. */
    size_t get_num_azimuths_per_elevation(size_t elevation) const {
        debug::assert2(elevation < m_num_elevations);
        return m_num_azimuths_per_elevation[elevation];
    }

    /** @brief Getter for num_elevations. */
    size_t get_num_elevations() const {
        return m_num_elevations;
    }

    /** @brief Getter for m_max_frequency. */
    float get_max_frequency() const {
        return m_max_frequency;
    }

    /** @brief Getter for m_mult_factor. */
    float get_mult_factor() const {
        return m_mult_factor;
    }

    /** @brief Getter for m_frequency_ratio. */
    float get_frequency_ratio() const {
        return m_frequency_ratio;
    }

    /** @brief Getter for m_angular_ratio. */
    float get_angular_ratio() const {
        return m_angular_ratio;
    }

    /** @brief Getter for m_lowpass_order. */
    float get_lowpass_order() const {
        return m_lowpass_order;
    }

    /** @brief Getter for m_lowpass_cutoff. */
    float get_lowpass_cutoff() const {
        return m_lowpass_cutoff;
    }

    /** @brief Getter for m_uniform_sampling. */
    bool get_uniform_sampling() const {
        return m_uniform_sampling;
    }

private:
    /**
     * @brief Creates a single log-Gabor filter.
     * 
     * @param[in] freq0  Central frequency value.
     * @param[in] phi0   Azimuth angle (in radians).
     * @param[in] theta0 Elevation angle (in radians).
     *
     * @returns The filter data array.
     */
    float* create_filter(float freq0,
                         float phi0,
                         float theta0);

    /**
     * @brief Reads a single log-Gabor filter from file.
     * 
     * @param[in] scale     Filter scale parameter.
     * @param[in] azimuth   Filter azimuth parameter.
     * @param[in] elevation Filter elevation parameter.
     *
     * @returns The filter data array.
     */
    float* read_filter(size_t scale,
                       size_t azimuth = 0,
                       size_t elevation = 0);

    /**
     * @brief Writes a single log-Gabor filter into a file.
     * 
     * @param[in] filter The filter data array.
     * @param[in] scale     Filter scale parameter.
     * @param[in] azimuth   Filter azimuth parameter.
     * @param[in] elevation Filter elevation parameter.
     */
    void write_filter(float *filter,
                      size_t scale,
                      size_t azimuth = 0,
                      size_t elevation = 0);

    /** @brief Setter for m_filename_prefix. */
    void set_filename_prefix(std::string filename_prefix) {
        debug::assert2(!filename_prefix.empty());
        m_filename_prefix = filename_prefix;
    }

    /** @brief Setter for m_sizes. */
    void set_sizes(triple<size_t> sizes) {
        m_sizes = sizes;
    }

    /** @brief Setter for m_num_scales. */
    void set_num_scales(size_t num_scales) {
        debug::assert2(num_scales > 0);
        m_num_scales = num_scales;
    }

    /** @brief Setter for m_num_azimuths. */
    void set_num_azimuths(size_t num_azimuths) {
        debug::assert2(num_azimuths > 0);
        m_num_azimuths = num_azimuths;
    }

    /** @brief Setter for m_num_elevations. */
    void set_num_elevations(size_t num_elevations) {
        debug::assert2(num_elevations > 0);
        m_num_elevations = num_elevations;
    }

    /** @brief Setter for m_max_frequency. */
    void set_max_frequency(float max_frequency) {
        debug::assert2(max_frequency > 0.0 && max_frequency < 0.5);
        m_max_frequency = max_frequency;
    }

    /** @brief Setter for m_mult_factor. */
    void set_mult_factor(float mult_factor) {
        debug::assert2(mult_factor > 0.0);
        m_mult_factor = mult_factor;
    }

    /** @brief Setter for m_frequency_ratio. */
    void set_frequency_ratio(float frequency_ratio) {
        debug::assert2(frequency_ratio > 0.0);
        m_frequency_ratio = frequency_ratio;
    }

    /** @brief Setter for m_angular_ratio. */
    void set_angular_ratio(float angular_ratio) {
        debug::assert2(angular_ratio > 0.0);
        m_angular_ratio = angular_ratio;
    }

    /** @brief Setter for m_lowpass_order. */
    void set_lowpass_order(float lowpass_order) {
        debug::assert2(lowpass_order > 0.0);
        m_lowpass_order = lowpass_order;
    }

    /** @brief Setter for m_lowpass_cutoff. */
    void set_lowpass_cutoff(float lowpass_cutoff) {
        debug::assert2(lowpass_cutoff > 0.0 && lowpass_cutoff <= 0.5);
        m_lowpass_cutoff = lowpass_cutoff;
    }

    /** @brief Setter for m_uniform_sampling. */
    void set_uniform_sampling(bool uniform_sampling) {
        m_uniform_sampling = uniform_sampling;
    }

    /** @attention Purposely not implemented. */
    log_gabor_filter_bank(const log_gabor_filter_bank &other);

    /** @attention Purposely not implemented. */
    log_gabor_filter_bank& operator=(const log_gabor_filter_bank &rhs);

private:
    std::string    m_filename_prefix;
    triple<size_t> m_sizes;
    size_t         m_num_scales;
    size_t         m_num_azimuths;
    size_t        *m_num_azimuths_per_elevation;
    size_t         m_num_elevations;
    float          m_max_frequency;
    float          m_mult_factor;
    float          m_frequency_ratio;
    float          m_angular_ratio;
    float          m_lowpass_order;
    float          m_lowpass_cutoff;
    bool           m_uniform_sampling;
};

// ----------------------------------------------------------------------------

/**
 * @brief Output stream operator for bip::log_gabor_filter_bank.
 */
std::ostream& operator<<(std::ostream &os,
                         const log_gabor_filter_bank &lgfb);


}

#endif
