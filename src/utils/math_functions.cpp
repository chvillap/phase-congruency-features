/**
 * @file   math_functions.cpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup utils
 * @ingroup    utils
 *
 * @copyright Copyright (c) 2016 Carlos Henrique Villa Pinto
 * @license GPL v2.0
 */

#include "math_functions.h"


namespace bip
{


const float EPSILON = 1E-5;


bool
fp_equal(float n1, float n2, float maxerror)
{
    // Try a simple equality test first.
    if (n1 == n2)
        return true;
    // Then try a "fuzzy" comparison using the absolute error.
    return fabs(n1 - n2) < maxerror;
}


triple<float>
cart2sph(triple<float> cart)
{
    const float x     = cart[0];
    const float y     = cart[1];
    const float z     = cart[2];
    const float rho   = sqrt(sqr(x) + sqr(y) + sqr(z));
    const float phi   = atan2(y, x);
    const float theta = atan2(z, sqrt(sqr(x) + sqr(y)));

    triple<float> sph = {rho, phi, theta};
    return sph;
}


triple<float>
sph2cart(triple<float> sph)
{
    const float rho   = sph[0];
    const float phi   = sph[1];
    const float theta = sph[2];
    const float x     = rho * cos(theta) * cos(phi);
    const float y     = rho * cos(theta) * sin(phi);
    const float z     = rho * sin(theta);

    triple<float> cart = {x, y, z};
    return cart;
}


void
FFT(fftwf_complex *M, triple<size_t> sizes, bool backward)
{
    const size_t dimensions = 3;

    const char WISDOM_FILENAME_FORWARD[32]  = "wisdom_fftwf_forward.txt";
    const char WISDOM_FILENAME_BACKWARD[32] = "wisdom_fftwf_backward.txt";

    fftwf_plan plan;
    FILE *wisdom_file = NULL;

    // Open the correct wisdom file according to the desired transform.
    if (backward)
        wisdom_file = fopen(WISDOM_FILENAME_BACKWARD, "r");
    else
        wisdom_file = fopen(WISDOM_FILENAME_FORWARD, "r");

    // Import FFTW plan settings from wisdom file if possible.
    // This will save a lot of runtime in this FFT computation.
    if (wisdom_file) {
        fftwf_import_wisdom_from_file(wisdom_file);
        fclose(wisdom_file);
    }

    // Reverse the order of the sizes array (so the size in Z axis is the first
    // element).
    int ft_sizes[3];
    for (size_t d = 0; d < dimensions; ++d)
        ft_sizes[d] = static_cast<int>(sizes[dimensions-d-1]);

    // Create the correct FFTW plan according to the desired transform.
    if (backward) {
        plan = fftwf_plan_dft(dimensions, ft_sizes, M, M, FFTW_BACKWARD,
                              FFTW_ESTIMATE);
    } else {
        plan = fftwf_plan_dft(dimensions, ft_sizes, M, M, FFTW_FORWARD,
                              FFTW_ESTIMATE);
    }

    // Export the FFTW plan settings to the wisdom file.
    // This will save a lot of runtime in future FFT computations.
    if (wisdom_file) {
        fftwf_export_wisdom_to_file(wisdom_file);
        fclose(wisdom_file);
    }

    // Compute the transform and clean up.
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);

    // The normalization step is applied only in backward FFT.
    if (backward) {
        size_t total_size = 1;
        for (size_t d = 0; d < dimensions; ++d)
            total_size *= sizes[d];
        
        for (size_t i = 0; i < total_size; ++i) {
            M[i][0] /= total_size;
            M[i][1] /= total_size;
        }
    }
}


float
hanning(triple<size_t> xyz, triple<size_t> sizes)
{
    const float hx = (sizes[0] == 1) ?
                     1.0 :
                     0.5 * (1.0 - cos(2*M_PI * xyz[0] / (sizes[0]-1)));
    const float hy = (sizes[1] == 1) ?
                     1.0 :
                     0.5 * (1.0 - cos(2*M_PI * xyz[1] / (sizes[1]-1)));
    const float hz = (sizes[2] == 1) ?
                     1.0 :
                     0.5 * (1.0 - cos(2*M_PI * xyz[2] / (sizes[2]-1)));
    return hx * hy * hz;
}


void
eigen(double *M, double *eigenvalues, double *eigenvectors)
{
    const size_t dimensions = 3;

    // Allocate GSL data.
    gsl_matrix_view m = gsl_matrix_view_array(M, dimensions, dimensions);
    gsl_vector *evals = gsl_vector_alloc(dimensions);
    gsl_matrix *evecs = gsl_matrix_alloc(dimensions, dimensions);
    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(dimensions);

    // Calculate the eigenvalues and eigenvectors of M.
    gsl_eigen_symmv(&m.matrix, evals, evecs, w);
    gsl_eigen_symmv_free(w);

    // Sort the eigenvectors in descending order by (signed) eigenvalues.
    gsl_eigen_symmv_sort(evals, evecs, GSL_EIGEN_SORT_VAL_DESC);

    // Set each eigenvector as a column of the resulting matrix.
    for (size_t i = 0; i < dimensions; ++i) {
        eigenvalues[i] = gsl_vector_get(evals, i);
        for (size_t j = 0; j < dimensions; ++j)
            eigenvectors[i*dimensions + j] = gsl_matrix_get(evecs, j, i);
    }
}


float
median(float *array, size_t size)
{
    debug::assert2(array != NULL);

    // Copy the data to an auxiliary vector.
    std::vector<float> aux;
    aux.assign(array, array + size);

    // Find the element in the middle of the data array.
    // Notice that this function rearranges the vector. That's why we need to
    // copy the original array first.
    size_t n = aux.size() / 2;
    std::nth_element(aux.begin(), aux.begin() + n, aux.end());

    // If the data size is odd, we have the median element.
    if(aux.size() % 2)
        return aux[n];
    
    // Otherwise, we compute the average of the nth and (n-1)th elements.
    std::nth_element(aux.begin(), aux.begin() + n-1, aux.end());
    return 0.5 * (aux[n] + aux[n-1]);
}


void
normalize_min_max(float *array, size_t size, float new_min, float new_max)
{
    float old_min =  FLT_MAX;
    float old_max = -FLT_MAX;

    // Find minimum and maximum.
    for (size_t i = 0; i < size; ++i) {
        if (old_min > array[i])
            old_min = array[i];
        if (old_max < array[i])
            old_max = array[i];
    }

    float old_range = old_max - old_min;
    float new_range = new_max - new_min;

    // Put data in the range [new_min, new_max].
    for (size_t i = 0; i < size; ++i)
        array[i] = (array[i] - old_min) / old_range;
    for (size_t i = 0; i < size; ++i)
        array[i] = array[i] * new_range + new_min;
}


}
