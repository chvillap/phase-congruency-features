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
     * Write it!
     */
}


// ------------------------------------------------------------------------------------------------

std::ostream&
operator<<(std::ostream &os, const phase_congruency &pc)
{
    return pc.print(os);
}


}
