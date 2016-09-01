/**
 * @file   phase_congruency.cpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup phase_congruency
 * @ingroup    phase_congruency
 *
 * @copyright Copyright (c) 2016 Carlos Henrique Villa Pinto
 * @license GPL v2.0
 */

#include "phase_congruency.h"


namespace bip
{


phase_congruency::
phase_congruency(std::string            filename_prefix,
                 float                 *input_image,
                 log_gabor_filter_bank *filter_bank,
                 triple<size_t>         sizes,
                 bool                  *input_mask,
                 float                  noise_threshold,
                 float                  noise_std,
                 float                  sigmoid_gain,
                 float                  sigmoid_cutoff)
{
    set_filename_prefix(filename_prefix);
    set_input_image(input_image);
    set_filter_bank(filter_bank);
    set_sizes(sizes);
    set_input_mask(input_mask);
    set_noise_threshold(noise_threshold);
    set_noise_std(noise_std);
    set_sigmoid_gain(sigmoid_gain);
    set_sigmoid_cutoff(sigmoid_cutoff);
}


phase_congruency::
~phase_congruency()
{
    // Nothing.
}


std::ostream&
phase_congruency::
print(std::ostream &os) const
{
    os << "{" 
       << "m_filename_prefix: " << m_filename_prefix << ", "
       << "m_filter_bank: "     << m_filter_bank     << ", "
       << "m_sizes: "           << m_sizes           << ", "
       << "m_input_image: "     << m_input_image     << ", "
       << "m_input_mask: "      << m_input_mask      << ", "
       << "m_noise_threshold: " << m_noise_threshold << ", "
       << "m_noise_std: "       << m_noise_std       << ", "
       << "m_sigmoid_gain: "    << m_sigmoid_gain    << ", "
       << "m_sigmoid_cutoff: "  << m_sigmoid_cutoff  << ", "
       << "}";

    return os;
}


void
phase_congruency::
compute()
{
    /*
     * TODO:
     * Parallelize this method by dividing the computations below in multiple
     * threads. There are many places this can be done. For example, each
     * log-Gabor filter is applied to the input image independently. Moreover,
     * many computations are done pixelwise, so the pixel access could be
     * parallelized too. The performance could be greatly improved with that.
     */

    // Get the total sizes in space and frequency domains.
    const size_t total_size   = m_sizes[0] * m_sizes[1] * m_sizes[2];
    const size_t f_total_size = total_size * sizeof(fftwf_complex);

    // Allocate data for the Fourier transform of the input image.
    fftwf_complex *f_input_image = (fftwf_complex*) fftwf_malloc(f_total_size);

    // Array of filtered images in the frequency domain.
    // It will contain one element per scale of the bank of filters.
    fftwf_complex *f_filtered_images = (fftwf_complex*)
        fftwf_malloc(f_total_size * m_filter_bank->get_num_scales());

    if (!f_input_image || !f_filtered_images)
        throw std::bad_alloc();

    // Compute the input image's Fourier transform (shifted).
    compute_shifted_FFT(f_input_image);


    // These arrays are used in the PC computation.
    float *sum_amplitude       = new float[total_size]();
    float *max_amplitude       = new float[total_size]();
    float *total_sum_amplitude = new float[total_size]();
    float *total_sum_energy    = new float[total_size]();
    float *cov_xx              = new float[total_size]();
    float *cov_xy              = new float[total_size]();
    float *cov_xz              = new float[total_size]();
    float *cov_yy              = new float[total_size]();
    float *cov_yz              = new float[total_size]();
    float *cov_zz              = new float[total_size]();
    float *pc_map              = new float[total_size]();
    float *directional_pc_map  = new float[total_size]();

    // And these are used to get the image features.
    float *moments_eigenvalues_maps[3] = {
        new float[total_size](),
        new float[total_size](),
        new float[total_size]()
    };
    triple<float> *moments_eigenvectors_maps[3] = {
        new triple<float>[total_size](),
        new triple<float>[total_size](),
        new triple<float>[total_size]()
    };
    triple<float> *directional_pc_max_map =
        new triple<float>[total_size]();
    

    #ifdef PHASE_CONGRUENCY_VERBOSE_ON
        std::cout << "Computing the PC maps...\n";
    #endif

    // Current orientation's index.
    size_t o = 0;

    // Get the step size in elevation angles.
    const float dtheta = (m_filter_bank->get_num_elevations() == 1) ?
                         0.0 :
                         M_PI_2 / (m_filter_bank->get_num_elevations() - 1);

    for (size_t e = 0; e < m_filter_bank->get_num_elevations(); ++e) {
        // Get the current elevation angle.
        const float theta     = e * dtheta;
        const float cos_theta = cos(theta);
        const float sin_theta = sin(theta);

        // Get the step size in azimuth angles.
        const float dphi = (m_filter_bank->get_num_azimuths() == 1) ?
                           0.0 :
                           (e == 0) ?
                               M_PI   / m_filter_bank->get_num_azimuths(0) :
                               M_PI*2 / m_filter_bank->get_num_azimuths(e);

        for (size_t a = 0; a < m_filter_bank->get_num_azimuths(e); ++a) {
            // Noise energy suppresion threshold.
            float T = 0.0;

            // Get the current azimuth angle.
            const float phi     = a * dphi;
            const float cos_phi = cos(phi);
            const float sin_phi = sin(phi);

            // The accumulated amplitudes are reset for each orientation.
            memset(sum_amplitude, 0, total_size * sizeof(float));
            memset(max_amplitude, 0, total_size * sizeof(float));

            for (size_t s = 0; s < m_filter_bank->get_num_scales(); ++s) {
                // Pointer to the filtered image at the current scale.
                fftwf_complex *p_f_filtered_image =
                    &f_filtered_images[s * total_size];

                // Apply a single log-Gabor filter (in the frequency domain) to
                // the input image. The result is stored in a slice of the
                // filtered images array.
                compute_filtering(p_f_filtered_image, f_input_image, s, a, e);

                // Accumulate amplitude responses over scales.
                for (size_t i = 0; i < total_size; ++i) {
                    // Ignore locations outside the region of interest.
                    if (m_input_mask && !m_input_mask[i])
                        continue;

                    const float even      = p_f_filtered_image[i][0];
                    const float odd       = p_f_filtered_image[i][1];
                    const float amplitude = sqrt(sqr(even) + sqr(odd));

                    sum_amplitude[i] += amplitude;
                    max_amplitude[i] = std::max(amplitude, max_amplitude[i]);
                }

                // Automatic noise energy threshold estimation.
                if (m_noise_threshold < 0.0 && s == 0) {
                    float tau = median(sum_amplitude, total_size) /
                                sqrt(log(4.0));
                    float invmult = 1.0 / m_filter_bank->get_mult_factor();
                    float nscales = m_filter_bank->get_num_scales();
                    float total_tau = tau * (1.0 - pow(invmult, nscales)) /
                                            (1.0 - invmult);
                    float noise_mean = total_tau * sqrt(M_PI_2);
                    float noise_sigma = total_tau * sqrt((4.0 - M_PI) / 2.0);

                    T = noise_mean + m_noise_std * noise_sigma;

                    #ifdef PHASE_CONGRUENCY_DEBUG_ON
                        std::cout << " (estimated T = " << T << ")";
                    #endif
                }

                #ifdef PHASE_CONGRUENCY_VERBOSE_ON
                    std::cout << " - done";
                #endif

                // The same block of memory is reused in the computation of all
                // directional PC maps, so data assigned in the previous
                // orientation must be cleaned.
                memset(directional_pc_map, 0, sizeof(directional_pc_map));

                /*
                 * TODO:
                 * Write the rest!
                 */
            }
        }
    }
}


void
phase_congruency::
compute_shifted_FFT(fftwf_complex *f_target)
{
    /*
     * TODO:
     * Parallelize this method by dividing the computation of the nested loop
     * below in multiple threads, since each pixel is computed independently of
     * the others. The performance could be greatly improved with that.
     */

    #ifdef PHASE_CONGRUENCY_VERBOSE_ON
        std::cout << "Computing the FFT";
    #endif

    // Allow the use of multithreading in the FFT computation.
    fftwf_init_threads();
    fftwf_plan_with_nthreads(DEFAULT_FFTW_NUMBER_OF_THREADS);

    // Shift the DC component to the center of the image.
    // The target is initialized with the result.
    for (size_t i = 0, z = 0; z < m_sizes[2]; ++z)
        for (size_t y = 0; y < m_sizes[1]; ++y)
            for (size_t x = 0; x < m_sizes[0]; ++x, ++i) {
                f_target[i][0] = m_input_image[i] * pow(-1.0, x + y + z);
                f_target[i][1] = 0.0;
            }

    // Compute the forward FFT of the input image.
    FFT(f_target, m_sizes);

    // Save the resulting frequency spectrum.
    #ifdef PHASE_CONGRUENCY_DEBUG_ON
    {
        size_t total_size = m_sizes[0] * m_sizes[1] * m_sizes[2];
        float *f_spectrum = new float[total_size]();

        for (size_t i = 0; i < total_size; ++i)
            f_spectrum[i] = log(1 + sqrt(sqr(f_target[i][0]) +
                                         sqr(f_target[i][1])));

        normalize_min_max(f_spectrum, total_size, 0.0, 1.0);

        /*
         * TODO:
         * Use some library to write 2D/3D images here.
         */
        delete[] f_spectrum;
    }
    #endif

    #ifdef PHASE_CONGRUENCY_VERBOSE_ON
        std::cout << " - done";
    #endif
}


void
phase_congruency::
compute_filtering(fftwf_complex *f_output,
                  fftwf_complex *f_input,
                  size_t         scale,
                  size_t         azimuth,
                  size_t         elevation)
{
    /*
     * TODO:
     * Parallelize this method by dividing the computation of the nested loops
     * below in multiple threads, since each pixel is computed independently of
     * the others. The performance could be greatly improved with that.
     */

    // Get a single log-Gabor filter (in the frequency domain) for
    // the given scale, azimuth and elevation.
    float *f_filter = m_filter_bank->get_filter(scale, azimuth, elevation);

    #ifdef PHASE_CONGRUENCY_VERBOSE_ON
        std::cout << "   Processing filter: "
                  << "sc = "
                  << std::setfill('0') << std::setw(3) << scale
                  << "az = "
                  << std::setfill('0') << std::setw(3) << azimuth
                  << "el = "
                  << std::setfill('0') << std::setw(3) << elevation;
    #endif

    // Apply the log-Gabor filter.
    for (size_t i = 0, z = 0; z < m_sizes[2]; ++z)
        for (size_t y = 0; y < m_sizes[1]; ++y)
            for (size_t x = 0; x < m_sizes[0]; ++x, ++i) {
                f_output[i][0] = f_filter[i] * f_input[i][0];
                f_output[i][1] = f_filter[i] * f_input[i][1];
            }

    delete[] f_filter;

    // Compute the backward FFT in order to get the filtered image in the
    // space domain.
    FFT(f_output, m_sizes, true);

    // Shift the DC component back to its original location.
    for (size_t i = 0, z = 0; z < m_sizes[2]; ++z)
        for (size_t y = 0; y < m_sizes[1]; ++y)
            for (size_t x = 0; x < m_sizes[0]; ++x, ++i) {
                f_output[i][0] *= pow(-1.0, x + y + z);
                f_output[i][1] *= pow(-1.0, x + y + z);
            }

    // Save the filter responses (even, odd, amplitude).
    #ifdef PHASE_CONGRUENCY_DEBUG_ON
    {
        const size_t total_size = m_sizes[0] * m_sizes[1] * m_sizes[2];

        float *f_even      = new float[total_size]();
        float *f_odd       = new float[total_size]();
        float *f_amplitude = new float[total_size]();

        for (size_t i = 0; i < total_size; ++i) {
            f_even[i]      = f_output[i][0];
            f_odd[i]       = f_output[i][1];
            f_amplitude[i] = sqrt(sqr(f_even[i]) + sqr(f_odd[i]));
        }
        normalize_min_max(f_even,      total_size, -1.0, 1.0);
        normalize_min_max(f_odd,       total_size, -1.0, 1.0);
        normalize_min_max(f_amplitude, total_size,  0.0, 1.0);

        char filename_suffix[16];

        sprintf(filename_suffix, "even_%02u_%02u_%02u.nii",
                scale, azimuth, elevation);
        /*
         * TODO:
         * Use some library to write 2D/3D images here.
         */

        sprintf(filename_suffix, "odd_%02u_%02u_%02u.nii",
                scale, azimuth, elevation);
        /*
         * TODO:
         * Use some library to write 2D/3D images here.
         */

        sprintf(filename_suffix, "amplitude_%02u_%02u_%02u.nii",
                scale, azimuth, elevation);
        /*
         * TODO:
         * Use some library to write 2D/3D images here.
         */

        delete[] f_even;
        delete[] f_odd;
        delete[] f_amplitude;
    }
    #endif
}


// ----------------------------------------------------------------------------

std::ostream&
operator<<(std::ostream &os, const phase_congruency &pc)
{
    return pc.print(os);
}


}
