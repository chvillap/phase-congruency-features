/**
 * @file   image_io.h
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup utils
 * @ingroup    utils
 *
 * @copyright Copyright (c) 2016 Carlos Henrique Villa Pinto
 * @license GPL v2.0
 */

#ifndef IMAGE_IO_H
#define IMAGE_IO_H

#include <cstddef>
#include <cstdlib>
#include <itkImage.h>
#include <itkVector.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include "assert2.h"
#include "types.h"


namespace bip
{
namespace io
{


/**
 * @fn read_image
 *
 * @brief Reads an image from file using ITK.
 *
 * @param[in] filename Image file name.
 *
 * @returns A smart pointer to the image object in ITK's format.
 */
template <class TPixel, size_t dims>
typename itk::Image<TPixel, dims>::Pointer read_image(std::string filename);


/**
 * @fn write_image
 *
 * @brief Writes an image using ITK.
 *
 * @param[in] filename Image file name.
 * @param[in] image    Image object in ITK's format.
 */
template <class TPixel, size_t dims>
void write_image(std::string filename,
                 typename itk::Image<TPixel, dims>::Pointer image);


/**
 * @fn array2image
 *
 * @brief Creates an ITK image from some array of pixel data.
 *
 * @param[in] data      Pixel data array.
 * @param[in] sizes     Number of array elements per dimension.
 * @param[in] reference Optional reference for 3D images.
 *
 * @returns A smart pointer to the image object in ITK's format.
 */
template <class TPixel, size_t dims>
typename itk::Image<TPixel, dims>::Pointer array2image(
    TPixel         *data,
    triple<size_t>  sizes,
    typename itk::Image<TPixel, dims>::Pointer reference = NULL);


/**
 * @fn image2array
 *
 * @brief Gets an array of pixel data from some ITK image.
 *
 * @param[in] image Image object in ITK's format.
 *
 * @returns A flattened array of pixel data.
 */
template <class TPixel, size_t dims>
TPixel* image2array(typename itk::Image<TPixel, dims>::Pointer image);


}
}

#include "image_io.cpp"

#endif
