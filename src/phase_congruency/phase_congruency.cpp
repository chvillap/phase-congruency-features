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
            // Noise energy threshold.
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

                // Accumulate amplitudes of filter responses over the scales.
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
                    T = estimate_noise_threshold(sum_amplitude);

                    #ifdef PHASE_CONGRUENCY_DEBUG_ON
                        std::cout << " (estimated T = " << T << ")";
                    #endif
                }

                #ifdef PHASE_CONGRUENCY_VERBOSE_ON
                    std::cout << " - done";
                #endif
            }

            // The same block of memory is reused in the computation of all
            // directional PC maps, so data assigned in the previous
            // orientation must be cleaned.
            memset(directional_pc_map, 0, sizeof(directional_pc_map));

            for (size_t i = 0; i < total_size; ++i) {
                // Ignore locations outside the region of interest.
                if (m_input_mask && !m_input_mask[i])
                    continue;

                // Accumulate the even and odd filter responses over scales.
                float sum_even = 0.0;
                float sum_odd  = 0.0;
                for (size_t s = 0; s < m_filter_bank->get_num_scales(); ++s)
                {
                    // Pointer to the filtered image at the current scale.
                    fftwf_complex *p_f_filtered_image =
                        &f_filtered_images[s * total_size];

                    // Accumulate the even and odd filter responses.
                    sum_even += p_f_filtered_image[i][0];
                    sum_odd  += p_f_filtered_image[i][1];
                }

                // Get the mean filter responses over scales.
                float norm      = sqrt(sqr(sum_even) + sqr(sum_odd));
                float mean_even = sum_even / (norm + EPSILON);
                float mean_odd  = sum_odd  / (norm + EPSILON);

                // Compute the local energy response for the current
                // orientation.
                float local_energy = 0.0;
                for (size_t s = 0; s < m_filter_bank->get_num_scales(); ++s)
                {
                    // Pointer to the filtered image at the current scale.
                    fftwf_complex *p_f_filtered_image =
                        &f_filtered_images[s * total_size];

                    float even = p_f_filtered_image[i][0];
                    float odd  = p_f_filtered_image[i][1];

                    // Theoretically, we need to compute the product of the
                    // amplitude of the filter responses by the phase deviation
                    // at the current pixel.
                    // In practice, this is the same as computing:
                    local_energy += (even * mean_even + odd * mean_odd) -
                                fabs(even * mean_odd - odd * mean_even);
                }

                // Apply the noise energy threshold to the local energy
                // (either an automatically calculated threshold or a manually
                // given one).
                if (m_noise_threshold < 0.0)
                    local_energy -= T;
                else
                    local_energy -= m_noise_threshold;

                // Apply the local energy weighting.
                local_energy = apply_energy_weighting(
                    local_energy, sum_amplitude[i], max_amplitude[i]);

                // Accumulate the total sums in amplitude and energy along
                // all orientations.
                total_sum_amplitude[i] += sum_amplitude[i];
                total_sum_energy[i]    += local_energy;

                // Compute the local phase congruency for the current location
                // and orientation.
                float local_pc = local_energy / (sum_amplitude[i] + EPSILON);

                // Set the pixel value for the current pixel of the current
                // orientation's directional PC map.
                directional_pc_map[i] = local_pc;

                // Update the pixel at the directional PC maxima map if needed.
                if (local_pc > directional_pc_max_map[i][0]) {
                    directional_pc_max_map[i][0] = local_pc;
                    directional_pc_max_map[i][1] = phi;
                    directional_pc_max_map[i][2] = theta;
                }

                // Project the local PC responses in the cartesian space.
                float proj_x = local_pc * cos_theta * cos_phi;
                float proj_y = local_pc * cos_theta * sin_phi;
                float proj_z = local_pc * sin_theta;

                // Accumulate the covariances between these projected responses
                // over the orientations.
                cov_xx[i] += sqr(proj_x);
                cov_yy[i] += sqr(proj_y);
                cov_zz[i] += sqr(proj_z);
                cov_xy[i] += proj_x * proj_y;
                cov_xz[i] += proj_x * proj_z;
                cov_yz[i] += proj_y * proj_z;
            }

            // Write the directional PC map.
            /*
             * TODO:
             * Use some library to write 2D/3D images here.
             */

            ++o; // Move on to the next orientation.
        }
    }

    // Write the directional PC maxima map.
    /*
     * TODO:
     * Use some library to write 2D/3D images here.
     */

    // Compute and write the final phase congruency map.
    for (size_t i = 0; i < total_size; ++i)
        pc_map[i] = total_sum_energy[i] / (total_sum_amplitude[i] + EPSILON);
    /*
     * TODO:
     * Use some library to write 2D/3D images here.
     */

    // Clean up the memory.
    delete[] sum_amplitude;
    delete[] max_amplitude;
    delete[] total_sum_amplitude;
    delete[] total_sum_energy;
    delete[] pc_map;
    delete[] directional_pc_map;
    delete[] directional_pc_max_map;
    fftwf_free(f_input_image);
    fftwf_free(f_filtered_images);
    fftwf_cleanup_threads();


    #ifdef PHASE_CONGRUENCY_VERBOSE_ON
        std::cout << "Computing PC moments, eigenvalues and eigenvectors";
    #endif

    // Covariance normalization factor.
    float orientations = static_cast<float>(
        m_filter_bank->get_num_orientations());
    float half_orientations = 0.5 * orientations;

    for (size_t i = 0, z = 0; z < m_sizes[2]; ++z)
        for (size_t y = 0; y < m_sizes[1]; ++y)
            for (size_t x = 0; x < m_sizes[0]; ++x, ++i) {
                // Ignore locations outside the region of interest.
                if (m_input_mask && !m_input_mask[i])
                    continue;

                // Finish the covariances calculation.
                cov_xx[i] /= half_orientations;
                cov_xy[i] /= orientations;
                cov_xz[i] /= orientations;
                cov_yy[i] /= half_orientations;
                cov_yz[i] /= orientations;
                cov_zz[i] /= half_orientations;

                // Compute the eigenvalues and eigenvectors of the covariance
                // matrix.
                // PS: if we have 1D or 2D data, all covariances related to Y
                // and/or Z are going to be 0, so it will be the same as
                // computing the eigenvalues/eigenvectors of a 1D or 2D matrix
                // in the end.
                double M[9] = {
                    cov_xx[i], cov_xy[i], cov_xz[i],
                    cov_xy[i], cov_yy[i], cov_yz[i],
                    cov_xz[i], cov_yz[i], cov_zz[i]
                };
                double eigenvalues[3];
                double eigenvectors[9];
                eigen(M, eigenvalues, eigenvectors);

                // Set the pixel at the eigenvalues and eigenvectors maps.
                for (size_t d = 0; d < 3; ++d) {
                    moments_eigenvalues_maps[d][i]     = eigenvalues[d];
                    moments_eigenvectors_maps[d][i][0] = eigenvectors[3*d];
                    moments_eigenvectors_maps[d][i][1] = eigenvectors[3*d + 1];
                    moments_eigenvectors_maps[d][i][2] = eigenvectors[3*d + 2];

                    // We are actually saving these maps as the eigenvectors
                    // scaled by their respective eigenvalues.
                    moments_eigenvectors_maps[d][i][0] *= eigenvalues[d];
                    moments_eigenvectors_maps[d][i][1] *= eigenvalues[d];
                    moments_eigenvectors_maps[d][i][2] *= eigenvalues[d];                    
                }
            }

    // Write all the eigenvalues and eigenvectors maps.
    for (size_t d = 0; d < 3; ++d) {
        char filename_suffix[16];

        sprintf(filename_suffix, "eigenvalues_%u", d);
        /*
        * TODO:
        * Use some library to write 2D/3D images here.
        */

        sprintf(filename_suffix, "eigenvectors_%u", d);
        /*
         * TODO:
         * Use some library to write 2D/3D images here.
         */
    }
    
    // Clean up memory.
    delete[] cov_xx;
    delete[] cov_xy;
    delete[] cov_xz;
    delete[] cov_yy;
    delete[] cov_yz;
    delete[] cov_zz;
    for (size_t d = 0; d < 3; ++d) {
        delete[] moments_eigenvalues_maps[d];
        delete[] moments_eigenvectors_maps[d];
    }

    #ifdef PHASE_CONGRUENCY_VERBOSE_ON
        std::cout << " - done";
    #endif
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


float
phase_congruency::
estimate_noise_threshold(float *sum_amplitude)
{
    const size_t total_size = m_sizes[0] * m_sizes[1] * m_sizes[2];

    float tau         = median(sum_amplitude, total_size) / sqrt(log(4.0));
    float invmult     = 1.0 / m_filter_bank->get_mult_factor();
    float nscales     = m_filter_bank->get_num_scales();
    float total_tau   = tau * (1.0 - pow(invmult, nscales)) / (1.0 - invmult);
    float noise_mean  = total_tau * sqrt(M_PI_2);
    float noise_sigma = total_tau * sqrt((4.0 - M_PI) / 2.0);

    return noise_mean + m_noise_std * noise_sigma;
}


float
phase_congruency::
apply_energy_weighting(float energy, float sum_amplitude, float max_amplitude)
{
    if (energy > 0.0) {
        // Get the frequency range width.
        // If there is only one non-zero component, width is 0.
        // If all components are equal, width is 1.
        float width = (sum_amplitude /  (max_amplitude + EPSILON) - 1.0) /
                        (m_filter_bank->get_num_scales() - 1);

        // The weighting function is a sigmoid.
        float weight = 1.0 + exp(m_sigmoid_gain * (m_sigmoid_cutoff - width));
        energy /= weight;
    }
    else // Negative weights are simply set to 0.
        energy = 0.0;

    return energy;
}


// ----------------------------------------------------------------------------

std::ostream&
operator<<(std::ostream &os, const phase_congruency &pc)
{
    return pc.print(os);
}


}
