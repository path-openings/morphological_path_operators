/*
 *		File:		ImageMagickIO.h
 *
 *		Purpose:	Wrap the ImageMagick file I/O functions
 *
 *		Author:		Ben Appleton

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

#ifndef IMAGEMAGICKIO_H
#define IMAGEMAGICKIO_H

#include <MagickCore/MagickCore.h>
#include "pde_toolbox_bimage.h"

/* Use dynamic (data-specific) constrast normalisation */
/* #define NORMALISE_CONTRAST_DYNAMIC */
/* Use known contrast range for our MRI data */
#define NORMALISE_CONTRAST_STATIC
#define MIN_MRI_PIXEL_VALUE 32768
#define MAX_MRI_PIXEL_VALUE 36864

/* Compiled against 16-bit version ImageMagick */
#define MIN_PIXEL_VALUE 0
#define MAX_PIXEL_VALUE 255
// Uncomment if you are using a 16-bit version of ImageMagick
// #define MAX_PIXEL_VALUE 65535

/*** Function prototypes ***/
/* Read a grayscale image and construct a BIMAGE object */
BIMAGE * read_grayscale_image(
    const char *currentpath,
	const char * filename
);

/* Write a BIMAGE object to a grayscale image */
void write_grayscale_image(
	BIMAGE * image,
	const char * filename
);

/* write_colour_image:
	Write a BIMAGE object to file as a colour image
*/
void write_colour_image(
	BIMAGE * bimage_red,
	BIMAGE * bimage_green,
	BIMAGE * bimage_blue,
	const char * filename
);

/* Normalise the image contrast */
void normalise_contrast(
	BIMAGE * image
);


#endif
