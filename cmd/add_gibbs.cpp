#include "fft.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <numeric>
#include <cmath>


#include "command.h"
#include "image.h"
#include "header.h"


using namespace MR;
using namespace App;

void usage() {
  AUTHOR = "Thea Bautista";

  SYNOPSIS = "Adding Gibbs-ringing to 3D image.";

  DESCRIPTION
    + "This reads an input image and outputs an image with added Gibbs-ringing.";

  ARGUMENTS
    + Argument ("inImg", "input image to be read").type_image_in()
    + Argument ("scale", "scale factor by which to downsample").type_integer(1, 20)
    + Argument ("outImg", "outuput image").type_image_out();

}



using ImageType = Image<cdouble>;



// gives proper index according to fourier indices
inline ssize_t indexshift(ssize_t n, ssize_t size) {
  if (n > size/2) n -= size;
  return n;
}



inline void set_index (const ImageType& lr, ImageType& hr, size_t axis)
{
  ssize_t pos = indexshift (lr.index(axis), lr.size(axis));
  hr.index(axis) = pos < 0 ? hr.size(axis) + pos : pos;
}


inline cdouble phase_shift (const ImageType& lr, const size_t factor)
{
  cdouble p (1.0, 0.0);
  for (int a = 0; a < 3; ++a)
    p *= std::exp (cdouble (0.0, ((factor-1.0)/factor) * (Math::pi*indexshift(lr.index(a), lr.size(a))) / lr.size(a)));
  return p;
}


void run()
{
  // reading input
  auto input = ImageType::open(argument[0]);

  // modifying header to create low resolution of input
  Header header (input);
  header.datatype() = DataType::CFloat32;
  auto highres_FT = ImageType::scratch (header, "FFT of input image");

  const int factor = argument[1];

  // create output image
  header.size(0) /= factor;
  header.size(1) /= factor;
  header.size(2) /= factor;
  header.spacing(0) *= factor;
  header.spacing(1) *= factor;
  header.spacing(2) *= factor;

  for (size_t j = 0; j < 3; ++j)
    for (size_t i = 0; i < 3; ++i)
      header.transform()(i,3) += 0.5 * (header.spacing(j) - input.spacing(j)) * header.transform()(i,j);

  auto lowres_FT = ImageType::scratch (header, "FFT of reduced image");

  header.datatype() = DataType::Float32;
  auto output = ImageType::create(argument[2],header);

  // 3D fft of input
  Math::FFT (input, highres_FT, 0, FFTW_FORWARD);
  Math::FFT (highres_FT, 1, FFTW_FORWARD);
  Math::FFT (highres_FT, 2, FFTW_FORWARD);

  for (auto l = Loop (lowres_FT, 0, 3) (lowres_FT); l; ++l) {
    set_index (lowres_FT, highres_FT, 0);
    set_index (lowres_FT, highres_FT, 1);
    set_index (lowres_FT, highres_FT, 2);
    lowres_FT.value() = phase_shift (lowres_FT, factor) * cdouble (highres_FT.value());
  }

  Math::FFT (lowres_FT, 0, FFTW_BACKWARD);
  Math::FFT (lowres_FT, 1, FFTW_BACKWARD);
  Math::FFT (lowres_FT, 2, FFTW_BACKWARD);

  const double scale = 1.0 / ( highres_FT.size(0) * highres_FT.size(1) * highres_FT.size(2));
  for (auto l = Loop (lowres_FT, 0, 3) (output, lowres_FT); l; ++l)
    output.value() = cdouble (lowres_FT.value()).real() * scale;
}


