/**
 * @file   image_io.cpp
 * @author Carlos H. Villa Pinto (chvillap@gmail.com)
 *
 * @addtogroup utils
 * @ingroup    utils
 *
 * @copyright Copyright (c) 2016 Carlos Henrique Villa Pinto
 * @license GPL v2.0
 */

namespace bip
{
namespace io
{


template <class TPixel, size_t dims>
typename itk::Image<TPixel, dims>::Pointer read_image(std::string filename)
{
    // Filename can't be empty.
    debug::assert2(!filename.empty());

    typedef itk::Image<TPixel, dims>     TImage;
    typedef itk::ImageFileReader<TImage> TReader;

    // Initialize the reader.
    typename TReader::Pointer reader = TReader::New();
    reader->SetFileName(filename.c_str());
    reader->Update();

    // Reads the file.
    return reader->GetOutput();
}


template <class TPixel, size_t dims>
void write_image(std::string filename,
                 typename itk::Image<TPixel, dims>::Pointer image)
{
    // Filename can't be empty and image can't be null.
    debug::assert2(!filename.empty());
    debug::assert2(image.IsNotNull());

    typedef itk::Image<TPixel, dims>     TImage;
    typedef itk::ImageFileWriter<TImage> TWriter;

    // Initialize the writer.
    typename TWriter::Pointer writer = TWriter::New();
    writer->SetInput(image);
    writer->SetFileName(filename.c_str());

    // Writes the file.
    writer->Update();
}


template <class TPixel, size_t dims>
typename itk::Image<TPixel, dims>::Pointer array2image(
    TPixel         *data,
    triple<size_t>  sizes,
    typename itk::Image<TPixel, dims>::Pointer reference)
{
    // Data array can't be null.
    debug::assert2(data != NULL);

    typedef itk::Image<TPixel, dims> TImage;

    // Set the image sizes in ITK format.
    typename TImage::SizeType sizes2;
    for (size_t d = 0; d < dims; ++d)
        sizes2[d] = sizes[d];

    // Set the image index.
    typename TImage::IndexType start;
    for (size_t d = 0; d < dims; ++d)
        start[d] = 0;

    // Set the image region.
    typename TImage::RegionType region;
    region.SetSize(sizes2);
    region.SetIndex(start);

    typename TImage::Pointer image = TImage::New();
    image->SetRegions(region);

    // Origin, spacing and direction can be taken from the reference image.
    if (reference.IsNotNull()) {
        image->SetOrigin(reference->GetOrigin());
        image->SetSpacing(reference->GetSpacing());
        image->SetDirection(reference->GetDirection());
    }

    // Initialize the pixel data.
    image->Allocate();
    image->FillBuffer(TPixel());

    // Set the image pixels.
    itk::ImageRegionIterator<TImage> iterator(image, image->GetBufferedRegion());
    size_t i = 0;
    while (!iterator.IsAtEnd()) {
        iterator.Set(data[i]);
        ++iterator;
        ++i;
    }

    return image;
}


template <class TPixel, size_t dims>
TPixel* image2array(typename itk::Image<TPixel, dims>::Pointer image)
{
    // Image can't be null.
    debug::assert2(image.IsNotNull());

    typedef itk::Image<TPixel, dims> TImage;

    // Get the image size.
    typename TImage::SizeType sizes = image->GetBufferedRegion().GetSize();
    size_t total_size = 1;
    for (size_t d = 0; d < dims; ++d)
        total_size *= sizes[d];

    // Allocate a flattened array of pixel data.
    TPixel *data = new TPixel[total_size]();

    // Get the image pixels.
    itk::ImageRegionConstIterator<TImage> iterator(image, image->GetBufferedRegion());
    size_t i = 0;
    while (!iterator.IsAtEnd()) {
        data[i] = iterator.Get();
        ++iterator;
        ++i;
    }

    return data;
}


}
}
