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


double inflection_point = 0.5;
bool elliptical = false;



void usage() {
  AUTHOR = "J-Donald Tournier (jdtournier@gmail.com)";

  SYNOPSIS = "Filter images with (1+cos(pi*k/kmax)/2 envelope in k-space";

  DESCRIPTION
    + "This reads an input nifti file and outputs an image after running fft function.";

  ARGUMENTS
    + Argument ("inImg", "input image to be read").type_image_in()
    + Argument ("outImg", "output image").type_image_out();

  OPTIONS
    + Option ("inflection_point", "set inflection point for Tukey window, "
        "as a floating-point value between 0 & 1, where 0 corresponds to a Hahn filter "
       " 1 to a square filter (default = "+str(inflection_point)+")")
    +  Argument ("value").type_float (0.0, 1.0)

    + Option ("elliptical", "also apply an elliptical filter, setting all values "
        "to zero outside of the largest ellipsoid that can fit in the acquired k-space");
}



using ImageType = Image<cdouble>;



// gives proper index according to fourier indices
inline double indexshift(ssize_t n, ssize_t size) {
  if (n > size/2) n -= size;
  return n;
}





// Filter for axis 0
class Filter {
  public:
    void operator() (ImageType& in){
      //apply filter
      cdouble val = in.value();
      double dist = 0.0;
      for (int n = 0; n < 3; ++n) {
        double x = 2.0 * std::abs (indexshift (in.index(n), in.size(n)) / in.size(n));
        dist += Math::pow2 (x);
        if (x > inflection_point)
          val *= 0.5 * (1.0 + std::cos (Math::pi * (x-inflection_point) / (1.0-inflection_point)));
      }
      if (elliptical && dist > 1.0) {
        in.value() = 0.0;
      }
      else {
        val /= in.size(0)*in.size(1)*in.size(2);
        in.value() = val;
      }
    }
};







void run()
{
  inflection_point = get_option_value ("inflection_point", inflection_point);
  elliptical = get_options ("elliptical").size() ? true : false;

  // reading input and assigning output
  auto input = ImageType::open(argument[0]);
  Header header (input);
  header.datatype() = DataType::CFloat32;
  auto image_FT = ImageType::scratch (header, "FFT of input image");

  auto output = ImageType::create (argument[1], header);


  { // full 3D FFT of input:
    ProgressBar progress ("performing Fourier transform", 3);
    Math::FFT (input, image_FT, 0, FFTW_FORWARD);
    ++progress;
    Math::FFT (image_FT, 1, FFTW_FORWARD);
    ++progress;
    Math::FFT (image_FT, 2, FFTW_FORWARD);
    ++progress;
  }

  ThreadedLoop("filtering in k-space", image_FT).run (Filter(), image_FT);

  { // full inverse 3D FFT of input:
    ProgressBar progress ("performing inverse Fourier transform", 3);
    Math::FFT (image_FT, 0, FFTW_BACKWARD);
    ++progress;
    Math::FFT (image_FT, 1, FFTW_BACKWARD);
    ++progress;
    Math::FFT (image_FT, output, 2, FFTW_BACKWARD);
    ++progress;
  }

}

