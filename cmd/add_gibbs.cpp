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
    + "This reads an input nifti file and outputs an image with added Gibbs-ringing.";

  ARGUMENTS
    + Argument ("inImg", "input image to be read").type_image_in()
    + Argument ("outImg", "outuput image").type_image_out();

}



using ImageType = Image<cdouble>;



// gives proper index according to fourier indices
inline double indexshift(ssize_t n, ssize_t size) {
  if (n > size/2) n -= size;
  return n;
}


void AddGibbsRinging (size_t axis, ImageType& input, ImageType& output) 
{
	int N = input.size(axis);

	for (int i = 0; i < N/8; i++) {
		output.index(axis) = input.index(axis) = i;
		output.value() = input.value();

		output.index(axis) = input.index(axis) = i + N;
		output.value() = input.value();
	}

}


void run()
{
	// reading input
	auto input = ImageType::open(argument[0]);
	
	// modifying header to create low resolution of input
	Header header (input);
	header.datatype() = DataType::CFloat32;
	auto image_FT = ImageType::scratch (header, "FFT of input image");
	// header.scale_() = 4.0;

	// create output image
	auto output = ImageType::create(argument[1],header);
	//auto image_FT = ImageType::scratch (header, "FFT of input image");

	// 3D fft of input
	Math::FFT(input,image_FT,0,FFTW_FORWARD);
	Math::FFT(image_FT,1,FFTW_FORWARD);
	Math::FFT(image_FT,2,FFTW_FORWARD);

	// add gibbs-ringing
	AddGibbsRinging (0, image_FT, output);
	AddGibbsRinging (1, image_FT, output);
	AddGibbsRinging (2, image_FT, output);

	Math::FFT (output, 0, FFTW_BACKWARD);
  Math::FFT (output, 1, FFTW_BACKWARD);
  Math::FFT (output, 2, FFTW_BACKWARD);
}


