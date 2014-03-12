/*
 *		File:		test_pathopen.cxx
 *
 *		Purpose:	Testing the pathopen implementation from Voir before modifying it
 *
 *		Author:		Ben Appleton
 *
 *              Changes :       Hugues Talbot 
 *
  Copyright Benjamin Appleton and Hugues Talbot, Nov 2009

  ben.appleton@gmail.com
  hugues.talbot@gmail.com / hugues.talbot@univ-paris-est.fr

This software is a computer program whose purpose is to perform
  connected morphological operators with path structuring elements.

This software is governed by the CeCILL-B  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL-B
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL-B license and that you accept its terms.

*/

#include <iostream>
#include <cstdlib>
#include <ctime>

using namespace std;

extern "C" {
	#include "pde_toolbox_bimage.h"
	#include "pde_toolbox_defs.h"
//	#include "pde_toolbox_LSTB.h"
	#include "ImageMagickIO.h"
}

#include "pathopenclose.h"

int usage(const char *name)
{
    cerr << "Usage : "<< name <<" <input image> L K <output image>" << endl;
    cerr << "Where : <input image> is an 8-bit grey-level image in any format readable by ImageMagick" << endl;
    cerr << "        L is the length of the path " << endl;
    cerr << "        K is the number of admissible missing pixels " << endl;
    cerr << "        <output image> is an 8-bit image. The extention determines the format" << endl;

    return 0;
}

int readargs(int argc, char *argv[], char **input, int *L, int *K, char **output)
{
    int notOKarg = 0;
    
    if (argc < 5) {
        notOKarg = 1;
    } else {
        *input = argv[1];
        *L = atoi(argv[2]);
        *K = atoi(argv[3]);
        *output = argv[4];
    }
    
    
    return notOKarg;
} 

int main(int argc, char **argv)
{
    int   i, L, K;
    char *input, *output;
    clock_t start, stop;
    
    if (readargs(argc, argv, &input, &L, &K, &output) != 0) {
        usage(argv[0]);
    } else {
	
	/* Open an image from file */
	BIMAGE * input_bimage = read_grayscale_image(argv[0], input);
	// Allocate remaining images
	BIMAGE * output_bimage = BIMAGE_constructor(input_bimage->dim);
	int nx = input_bimage->dim->buf[0];
	int ny = input_bimage->dim->buf[1];
        int num_pixels = nx*ny;
	PATHOPEN_PIX_TYPE * input_image = new PATHOPEN_PIX_TYPE[nx * ny];
	PATHOPEN_PIX_TYPE * output_image = new PATHOPEN_PIX_TYPE[nx * ny];

	// Convert intermediate float to PATHOPEN_PIX_TYPE (unsigned char)
	for (i = 0; i < num_pixels; ++i) {
            input_image[i] = static_cast<PATHOPEN_PIX_TYPE>(input_bimage->buf[i]);
	}

	cout << "Calling pathopen()" << endl;
        start = clock();
	pathopen(
            input_image, /* The input image */
            nx, ny,	 /* Image dimensions */
            L,		 /* The threshold line length */
            K,		 /* The maximum number of gaps in the path */
            output_image /* Output image */
            );
        stop = clock();
	cout << "pathopen() returned! CPU time elapsed:" << ((double)stop-start)/CLOCKS_PER_SEC << endl;

	/* Save output to file */
	// Convert 0-65535 float to PATHOPEN_PIX_TYPE (unsigned char)
	for (i = 0; i < num_pixels; ++i) {
            output_bimage->buf[i] = static_cast<PATHOPEN_PIX_TYPE>(output_image[i]);
	}
	// Write file
	write_grayscale_image(
		output_bimage,
		output
	);

	// Deallocate
	delete input_image;
	delete output_image;
	BIMAGE_destructor(input_bimage);
	BIMAGE_destructor(output_bimage);
        
    }
        
	return 0;
}
