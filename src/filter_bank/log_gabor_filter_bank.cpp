/**
 * @file   log_gabor_filter_bank.cpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup filter_bank
 * @ingroup    filter_bank
 *
 * @copyright Copyright (c) 2016 Carlos Henrique Villa Pinto
 * @license GPL v2.0
 */

#include "log_gabor_filter_bank.h"


namespace bip
{


log_gabor_filter_bank::
log_gabor_filter_bank(std::string    filename_prefix,
                      triple<size_t> sizes,
                      size_t         num_scales,
                      size_t         num_azimuths,
                      size_t         num_elevations,
                      float          max_frequency,
                      float          mult_factor,
                      float          frequency_ratio,
                      float          angular_ratio,
                      float          lowpass_order,
                      float          lowpass_cutoff,
                      bool           uniform_sampling)
{
    set_filename_prefix(filename_prefix);
    set_sizes(sizes);
    set_num_scales(num_scales);
    set_num_azimuths(num_azimuths);
    set_num_elevations(num_elevations);
    set_max_frequency(max_frequency);
    set_mult_factor(mult_factor);
    set_frequency_ratio(frequency_ratio);
    set_angular_ratio(angular_ratio);
    set_lowpass_order(lowpass_order);
    set_lowpass_cutoff(lowpass_cutoff);
    set_uniform_sampling(uniform_sampling);

    m_num_azimuths_per_elevation = new size_t[m_num_elevations]();
    m_num_azimuths_per_elevation[0] = m_num_azimuths;

    // With uniform sampling, every elevation has the same number of azimuth
    // angles. With nonuniform sampling, we reduce the number of azimuths as
    // the elevation angle increases, according to the formula below.
    if (m_uniform_sampling) {
        for (size_t e = 1; e < m_num_elevations; ++e)
            m_num_azimuths_per_elevation[e] = m_num_azimuths;
    } else {
        if (m_num_elevations > 1)
            m_num_azimuths_per_elevation[m_num_elevations-1] = 1;

        for (size_t e = 1; e < m_num_elevations-1; ++e) {
            float n = m_num_azimuths * cos(e * M_PI_2 / (m_num_elevations-1));
            m_num_azimuths_per_elevation[e] = (size_t) round(n) * 2;
        }
    }
}


log_gabor_filter_bank::
~log_gabor_filter_bank()
{
    delete[] m_num_azimuths_per_elevation;
}


std::ostream&
log_gabor_filter_bank::
print(std::ostream &os) const
{
    os << "{"
       << "m_filename_prefix: "  << m_filename_prefix  << ", "
       << "m_sizes: "            << m_sizes            << ", "
       << "m_num_scales: "       << m_num_scales       << ", "
       << "m_num_azimuths: "     << m_num_azimuths     << ", "
       << "m_num_elevations: "   << m_num_elevations   << ", "
       << "m_max_frequency: "    << m_max_frequency    << ", "
       << "m_mult_factor: "      << m_mult_factor      << ", "
       << "m_frequency_ratio: "  << m_frequency_ratio  << ", "
       << "m_angular_ratio: "    << m_angular_ratio    << ", "
       << "m_lowpass_order: "    << m_lowpass_order    << ", "
       << "m_lowpass_cutoff: "   << m_lowpass_cutoff   << ", "
       << "m_uniform_sampling: " << (m_uniform_sampling ? "true" : "false")
       << "}";

    return os;
}


float*
log_gabor_filter_bank::
get_filter(size_t scale, size_t azimuth, size_t elevation)
{
    debug::assert2(elevation < m_num_elevations);
    debug::assert2(azimuth < m_num_azimuths_per_elevation[elevation]);
    debug::assert2(scale < m_num_scales);

    #ifdef LOG_GABOR_FILTER_BANK_DEBUG_ON
        std::cout << "Getting filter: "
                  << "sc = "
                  << std::setfill('0') << std::setw(3) << scale
                  << "az = "
                  << std::setfill('0') << std::setw(3) << azimuth
                  << "el = "
                  << std::setfill('0') << std::setw(3) << elevation;
    #endif

    // Try to read an existing filter from file.
    float *filter = read_filter(scale, azimuth, elevation);

    #ifdef LOG_GABOR_FILTER_BANK_DEBUG_ON
        std::cout << " - done\n";
    #endif

    return filter;
}


void
log_gabor_filter_bank::
compute()
{
    /*
     * TODO:
     * Parallelize this method by dividing the computation of the nested loops
     * below in multiple threads. Each filter is computed independently of the
     * others, so the performance could be greatly improved with that.
     *
     * PS: choose between parallelizing this method or create_filter() method.
     * Parallelizing both is probably not a good idea (thread overhead).
     */

    // Variation in the elevation angle.
    float dtheta = (m_num_elevations == 1) ?
                   0.0 :
                   M_PI_2 / (m_num_elevations - 1);

    for (size_t e = 0; e < m_num_elevations; ++e) {
        // Variation in the azimuth angle.
        float dphi = 0.0;
        if (m_num_azimuths > 1) {
            dphi = (e == 0) ?
                   M_PI / m_num_azimuths_per_elevation[e] :
                   M_PI / m_num_azimuths_per_elevation[e] * 2;
        }
        // Central elevation angle.
        float theta0 = e * dtheta;

        for (size_t a = 0; a < m_num_azimuths_per_elevation[e]; ++a) {
            // Central azimuth angle.
            float phi0 = a * dphi;

            for (size_t s = 0; s < m_num_scales; ++s) {
                // Central frequency.
                float freq0 = m_max_frequency / pow(m_mult_factor, s);

                #ifdef LOG_GABOR_FILTER_BANK_VERBOSE_ON
                    std::cout << "Creating filter: "
                              << "sc = "
                              << std::setfill('0') << std::setw(3) << s
                              << ", az = "
                              << std::setfill('0') << std::setw(3) << a
                              << ", el = "
                              << std::setfill('0') << std::setw(3) << e;
                #endif

                // Create the filter for the current frequency and orientation.
                float *filter = create_filter(freq0, phi0, theta0);
                write_filter(filter, s, a, e);
                delete[] filter;

                #ifdef LOG_GABOR_FILTER_BANK_VERBOSE_ON
                    std::cout << " - done\n";
                #endif
            }
        }
    }
}


log_gabor_filter_bank*
log_gabor_filter_bank::
read_parameters(std::string filename)
{
    debug::assert2(!filename.empty());
    debug::assert2(filename.find_last_of(".bof") != std::string::npos);

    triple<size_t> sizes;
    size_t num_scales, num_azimuths, num_elevations;
    float max_frequency, mult_factor;
    float frequency_ratio, angular_ratio;
    float lowpass_order, lowpass_cutoff;
    bool  uniform_sampling;

    #ifdef LOG_GABOR_FILTER_BANK_VERBOSE_ON
        std::cout << "Reading parameters of the bank of filters from: "
                  << filename;
    #endif

    // Open a binary file input stream.
    std::fstream ifs(filename.c_str(), std::ios::in | std::ios::binary);
    if (ifs.fail())
        throw "read_parameters() : Couldn't open file stream";

    // Read the filter parameters.
    ifs.read(reinterpret_cast<char*>(&sizes[0]),         sizeof(size_t));
    ifs.read(reinterpret_cast<char*>(&sizes[1]),         sizeof(size_t));
    ifs.read(reinterpret_cast<char*>(&sizes[2]),         sizeof(size_t));
    ifs.read(reinterpret_cast<char*>(&num_scales),       sizeof(size_t));
    ifs.read(reinterpret_cast<char*>(&num_azimuths),     sizeof(size_t));
    ifs.read(reinterpret_cast<char*>(&num_elevations),   sizeof(size_t));
    ifs.read(reinterpret_cast<char*>(&max_frequency),    sizeof(float));
    ifs.read(reinterpret_cast<char*>(&mult_factor),      sizeof(float));
    ifs.read(reinterpret_cast<char*>(&frequency_ratio),  sizeof(float));
    ifs.read(reinterpret_cast<char*>(&angular_ratio),    sizeof(float));
    ifs.read(reinterpret_cast<char*>(&lowpass_order),    sizeof(float));
    ifs.read(reinterpret_cast<char*>(&lowpass_cutoff),   sizeof(float));
    ifs.read(reinterpret_cast<char*>(&uniform_sampling), sizeof(bool));
    ifs.close();

    // The filename prefix is everything that comes before "_#_".
    std::string filename_prefix(filename.substr(0, filename.find("_#_")));

    // Allocate the bank of filters.
    log_gabor_filter_bank *bof =
        new log_gabor_filter_bank(filename_prefix, sizes, num_scales,
                                  num_azimuths, num_elevations, max_frequency,
                                  mult_factor, frequency_ratio, angular_ratio,
                                  lowpass_order, lowpass_cutoff,
                                  uniform_sampling);

    #ifdef LOG_GABOR_FILTER_BANK_VERBOSE_ON
        std::cout << " - done\n";
    #endif

    return bof;
}


void
log_gabor_filter_bank::
write_parameters(log_gabor_filter_bank &bof)
{
    char filename_suffix[64];

    sprintf(filename_suffix,
            "_#_%03u_%03u_%03u_%02u_%02u_%02u_%03u_%03u_%03u_%03u_%1u.bof",
            bof.get_sizes()[0],
            bof.get_sizes()[1],
            bof.get_sizes()[2],
            bof.get_num_scales(),
            bof.get_num_azimuths(),
            bof.get_num_elevations(),
            static_cast<size_t>(bof.get_max_frequency()   * 100),
            static_cast<size_t>(bof.get_mult_factor()     * 100),
            static_cast<size_t>(bof.get_frequency_ratio() * 100),
            static_cast<size_t>(bof.get_angular_ratio()   * 100),
            static_cast<size_t>(bof.get_uniform_sampling()));

    std::string filename(bof.get_filename_prefix() +
                         std::string(filename_suffix));

    #ifdef LOG_GABOR_FILTER_BANK_VERBOSE_ON
        std::cout << "Writing parameters of the bank of filters to: "
                  << filename;
    #endif

    triple<size_t> sizes            = bof.get_sizes();
    size_t         scales           = bof.get_num_scales();
    size_t         azimuths         = bof.get_num_azimuths();
    size_t         elevations       = bof.get_num_elevations();
    float          max_frequency    = bof.get_max_frequency();
    float          mult_factor      = bof.get_mult_factor();
    float          frequency_ratio  = bof.get_frequency_ratio();
    float          angular_ratio    = bof.get_angular_ratio();
    float          lowpass_order    = bof.get_lowpass_order();
    float          lowpass_cutoff   = bof.get_lowpass_cutoff();
    bool           uniform_sampling = bof.get_uniform_sampling();

    // Open a binary file output stream.
    std::fstream ofs(filename.c_str(), std::ios::out | std::ios::binary);
    if (ofs.fail())
        throw "write_parameters() : Couldn't open file stream";

    // Write the filter parameters.
    ofs.write(reinterpret_cast<char*>(&sizes[0]),         sizeof(size_t));
    ofs.write(reinterpret_cast<char*>(&sizes[1]),         sizeof(size_t));
    ofs.write(reinterpret_cast<char*>(&sizes[2]),         sizeof(size_t));
    ofs.write(reinterpret_cast<char*>(&scales),           sizeof(size_t));
    ofs.write(reinterpret_cast<char*>(&azimuths),         sizeof(size_t));
    ofs.write(reinterpret_cast<char*>(&elevations),       sizeof(size_t));
    ofs.write(reinterpret_cast<char*>(&max_frequency),    sizeof(float));
    ofs.write(reinterpret_cast<char*>(&mult_factor),      sizeof(float));
    ofs.write(reinterpret_cast<char*>(&frequency_ratio),  sizeof(float));
    ofs.write(reinterpret_cast<char*>(&angular_ratio),    sizeof(float));
    ofs.write(reinterpret_cast<char*>(&lowpass_order),    sizeof(float));
    ofs.write(reinterpret_cast<char*>(&lowpass_cutoff),   sizeof(float));
    ofs.write(reinterpret_cast<char*>(&uniform_sampling), sizeof(bool));
    ofs.close();

    #ifdef LOG_GABOR_FILTER_BANK_VERBOSE_ON
        std::cout << " - done\n";
    #endif
}


float*
log_gabor_filter_bank::
create_filter(float freq0, float phi0, float theta0)
{
    /*
     * TODO:
     * Parallelize this method by dividing the computation of the nested loops
     * below in multiple threads. Each voxel is computed independently of the
     * others, so the performance could be greatly improved with that.
     */

    // Angular standard deviation.
    const float std_alpha = M_PI / m_num_azimuths / m_angular_ratio;

    // Trigonometric constants.
    const float sin_theta         = sin(theta0);
    const float cos_theta         = cos(theta0);
    const float cos_theta_cos_phi = cos_theta * cos(phi0);
    const float cos_theta_sin_phi = cos_theta * sin(phi0);

    // Allocate the filter data array.
    float *filter = new float[m_sizes[0] * m_sizes[1] * m_sizes[2]]();

    // Iterate through the frequency domain coordinates.
    for (size_t i = 0, z = 0; z < m_sizes[2]; ++z) {
        float w = (m_sizes[2] == 1) ?
                  0.0 :
                  0.5 - (float) z / m_sizes[2];

        for (size_t y = 0; y < m_sizes[1]; ++y) {
            float v = (m_sizes[1] == 1) ?
                      0.0 :
                      0.5 - (float) y / m_sizes[1];

            for (size_t x = 0; x < m_sizes[0]; ++x, ++i) {
                float u = (m_sizes[0] == 1) ?
                          0.0 :
                         -0.5 + (float) x / m_sizes[0];

                // Get the frequency value.
                float freq = sqrt(sqr(u) + sqr(v) + sqr(w));

                if (freq < EPSILON)
                    filter[i] = 0.0;
                else {
                    // Project the sample point in the cartesian space.
                    float uu = u * cos_theta_cos_phi;
                    float vv = v * cos_theta_sin_phi;
                    float ww = w * sin_theta;

                    // Angular distance between the central frequency point
                    // and the current sample point.
                    float alpha = acos((uu + vv + ww) / (freq + EPSILON));

                    // Angular and radial components in the frequency domain.
                    float spread = exp(-0.5 * sqr(alpha / std_alpha));
                    float radius = exp(-0.5 * sqr(log(freq / freq0) /
                                                  log(m_frequency_ratio)));

                    // // Angular distance between the central frequency point
                    // // and the current sample.
                    // float alpha = fabs(acos((uu + vv + ww) / (freq + EPSILON)));
                    // alpha = std::min(alpha * m_num_azimuths/2.0, M_PI);

                    // // Angular and radial components in the frequency domain.
                    // float spread = 0.5 * (1 + cos(alpha));
                    // float radius = exp(-0.5 * sqr(log(freq / freq0) /
                    //                               log(m_frequency_ratio)));

                    // Apply a Butterworth low-pass filter.
                    radius /= (1 + pow((freq / m_lowpass_cutoff),
                                       m_lowpass_order * 2));

                    // Combine the radial and angular spread components.
                    filter[i] = radius * spread;
                }
            }
        }
    }
    return filter;
}


float*
log_gabor_filter_bank::
read_filter(size_t scale, size_t azimuth, size_t elevation)
{
    char filename_suffix[64];

    // Set the full filename.
    sprintf(filename_suffix, "_#_%03u_%03u_%03u_%02u_%02u_%02u.nii",
            m_sizes[0], m_sizes[1], m_sizes[2], scale, azimuth, elevation);

    std::string filename(m_filename_prefix + std::string(filename_suffix));

    // Read the filter image from a file.
    float *filter = io::image2array<float, 3>(
                        io::read_image<float, 3>(filename));

    return filter;
}


void
log_gabor_filter_bank::
write_filter(float *filter, size_t scale, size_t azimuth, size_t elevation)
{
    char filename_suffix[64];

    // Set the full filename.
    sprintf(filename_suffix, "_#_%03u_%03u_%03u_%02u_%02u_%02u.nii",
            m_sizes[0], m_sizes[1], m_sizes[2], scale, azimuth, elevation);

    std::string filename(m_filename_prefix + std::string(filename_suffix));

    // Write the filter image into a file.
    io::write_image<float, 3>(
        filename, io::array2image<float, 3>(filter, m_sizes));
}


// ----------------------------------------------------------------------------

std::ostream&
operator<<(std::ostream &os, const log_gabor_filter_bank &lgfb)
{
    return lgfb.print(os);
}


}
