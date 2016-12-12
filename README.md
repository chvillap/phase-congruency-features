2D/3D phase congruency features
================
C++ implementation of the phase congruency technique for 2D and 3D data.

**Phase congruency** is a signal processing technique better known for its use in the detection of invariant image features. Rather than assuming that an image should be compressed into a set of edges, the phase congruency model assumes that the compressed image format should be high in information (or entropy), and low in redundancy. Thus, instead of searching for points where there are sharp changes in intensity, this model searches for patterns of order in the phase component of the Fourier transform. Phase is chosen because experiments demonstrated that it is crucial to the perception of visual features. Further physiological evidence indicates that the human visual system responds strongly to points in an image where the phase information is highly ordered. Thus the phase congruency model defines features as points in an image with high phase order. And for that it has a series of advantages over other image feature detectors.

The phase congruency measure is proportional to the local energy of the signal, therefore it can be calculated via convolution of the original image with a bank of spatial filters in quadrature. A bank of [log-Gabor filters](https://en.wikipedia.org/wiki/Log_Gabor_filter) is especially suited for that, thus an implementation of 2D/3D log-Gabor filters is also part of this project.

More details about phase congruency and its applications can be found in the following papers:

> Kovesi, P., 2000. Phase congruency: A low-level image invariant. Psychological Research 64, 136-148.

> Kovesi, P., 2003. Phase congruency detects corners and edges, in: The Australian Pattern Recognition Society Conference: DICTA, pp. 309-318.

> Ferrari, R.J., Allaire, S., Hope, A., Kim, J., Jaffray, D., Pekar, V., 2011. Detection of point landmarks in 3D medical images via phase congruency model. Journal of the Brazilian Computer Society 17, 117-132.

> Villa Pinto, C.H.; Ferrari, R.J., 2016. Initialization of deformable models in 3D magnetic resonance images guided by automatically detected phase congruency point landmarks. Pattern Recognition Letters 79, 1-7.

...among several others. In addition, [Dr. Peter Kovesi's website](http://www.peterkovesi.com) contains some great MATLAB implementations for 2D images.

## Dependencies

- [CMake 3.6](https://cmake.org)
- [FFTW 3.3](http://www.fftw.org)
- [GSL 2.2](https://www.gnu.org/software/gsl)
- [ITK 4.10](https://www.itk.org)

## Notes

- This implementation does not make use of [monogenic filters](https://www.math.ucdavis.edu/~saito/data/phase2/monogenic.pdf) due to the application it is aimed to. So be aware that this is certainly not the fastest implementation that can be achieved.
- The NIfTI (.nii) format is used for most image outputs. You can use softwares like [3D Slicer](https://www.slicer.org) and [ITK-Snap](http://www.itksnap.org/pmwiki/pmwiki.php) to open such files.
