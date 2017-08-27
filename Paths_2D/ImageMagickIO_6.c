/*
 *		File:		ImageMagickIO.c
 *
 *		Purpose:	Wrap the ImageMagick file I/O functions
 *
 *		Author:		Ben Appleton
 *
 *              Modified by:    Hugues Talbot	30 Nov 2009
 *                              Hugues Talbot   08 Aug 2017: Magick v7
 *
 *		Date:		05/05/2005


  Copyright Benjamin Appleton and Hugues Talbot, May 2005

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
 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <MagickCore/MagickCore.h>
#include "pde_toolbox_bimage.h"
#include "ImageMagickIO.h"

/* read_grayscale_image:
   Read a grayscale image from the specified file
   Construct and return a BIMAGE object
*/
BIMAGE *read_grayscale_image(const char *currentpath, const char * filename)
{
    BIMAGE * bimage;
    BVECT * bvect;
    int i, num_pixels;

    /* Various objects needed by ImageMagick */
    ExceptionInfo exception;
    Image * image;
    ImageInfo * image_info;
    PixelPacket * pixel_packet;

    /* Initialise ImageMagick.  This just tells it what directory (*argv) we're running from */
    //InitializeMagick((char *)NULL);
    MagickCoreGenesis(currentpath, MagickFalse);

    /* Initialise the ExceptionInfo - ImageMagick's way of implementing exceptions in C */
    GetExceptionInfo(&exception);

    /* Construct an image_info object by 'cloning' a null image_info object */
    image_info = CloneImageInfo((ImageInfo *) NULL);

    /*** Read an image from file ***/
    /* Set the name of the input file, which is stored in image_info */
    (void) strcpy(image_info->filename, filename);
    /* Actually read the image */
    image = ReadImage(image_info, &exception);
    /* Deal with exceptions (ie. drop out) */
    if (exception.severity != UndefinedException)
        CatchException(&exception);
    if (image == (Image *) NULL)
        exit(1);

    /*** Construct a BIMAGE ***/
    /* Setup the dimensions */
    bvect = BVECT_constructor(2);
    bvect->buf[0] = image->columns;
    bvect->buf[1] = image->rows;

    /* Construct an image of the desired dimensions */
    bimage = BIMAGE_constructor(bvect);
	
    /* Destroy the dimensions */
    BVECT_destructor(bvect);

    /*** Copy the image data into the BIMAGE ***/
    num_pixels = BVECT_prod(bimage->dim);
    pixel_packet = GetAuthenticPixels(image, 0, 0, image->columns, image->rows, &exception);
    if (pixel_packet == NULL) {
        printf("read_grayscale_image: GetImagePixels failed!\n");
        return (BIMAGE * )NULL;
    }
    for (i = 0; i < num_pixels; ++i) {
        bimage->buf[i] = (float)(pixel_packet[i].red + pixel_packet[i].green + pixel_packet[i].blue)/3.0f;
    }

    /*** Cleanup ***/
    /* Destroy the image info object */
    image_info = DestroyImageInfo(image_info);

    DestroyImage(image);

    /* Finalise the ExceptionInfo */
    DestroyExceptionInfo(&exception);
    
    /* Finalise ImageMagick */
    MagickCoreTerminus();
    
    //normalise_contrast(bimage);

    return bimage;
}


/* write_grayscale_image:
   Write a BIMAGE object to file as a grayscale image
*/
void write_grayscale_image(
    BIMAGE * bimage,
    const char * filename
    )
{
    /* Various objects needed by ImageMagick */
    ExceptionInfo exception;
    Image * image;
    ImageInfo * image_info;
    PixelPacket * pixel_packet;

    int i, num_pixels, maxval=0;

    /* Initialise ImageMagick.  This just tells it what directory (*argv) we're running from */
    MagickCoreGenesis((char *)NULL, MagickFalse);

    /* Initialise the ExceptionInfo - ImageMagick's way of implementing exceptions in C */
    GetExceptionInfo(&exception);

    /* Construct an image_info object by 'cloning' a null image_info object */
    image_info = CloneImageInfo((ImageInfo *) NULL);

    /*** Construct an ImageMagick Image of the appropriate dimensions ***/
    image = ConstituteImage(
        bimage->dim->buf[0],			// nx
        bimage->dim->buf[1],			// ny
        "I",							// Grayscale image, contains only "Intensity"
        CharPixel,						// Data type
        bimage->buf,					// Image data
        &exception						// Exception object
	);

    /*** Copy pixel values from BIMAGE object to ImageMagick's Image object ***/
    num_pixels = BVECT_prod(bimage->dim);
    pixel_packet = GetAuthenticPixels(image, 0, 0, image->columns, image->rows, &exception);
    if (pixel_packet == NULL) {
        printf("write_grayscale_image: GetImagePixels failed!\n");
        return;
    }

    /* greyscale -> repeated color channels */
    for (i = 0; i < num_pixels; ++i) {
        unsigned short int val = (unsigned short int)bimage->buf[i];
        if (val > maxval)
            maxval = val;
        pixel_packet[i].red = pixel_packet[i].green = pixel_packet[i].blue = val;
    }
    /* Synchronise the buffer with the contour image!! */
    SyncAuthenticPixels(image, &exception);

    /*** Write the image to file ***/
    /* To save an image to a specific location, store the filename in the image itself */
    (void) strcpy(image->filename, filename);
    /* Actually write the image */

    //image->depth=8; // force 8-bit output if possible
    if ((GetImageDepth(image,&exception) > 8) && (maxval <= 255)) {
            fprintf(stderr, "Reducing output image depth\n");
        for (i = 0 ; i < num_pixels ; ++i) {
            pixel_packet[i].red <<= 8;
            pixel_packet[i].green <<= 8;
            pixel_packet[i].blue <<= 8;
        }
        /* this actually truncates values */
        SetImageDepth(image, 8);
    }
    WriteImage(image_info, image);

    DestroyImage(image);

    /* Destroy the image info object */
    image_info = DestroyImageInfo(image_info);
    
    /* Finalise the ExceptionInfo */
    DestroyExceptionInfo(&exception);
    
    /* Finalise ImageMagick */
    MagickCoreTerminus();
    
    return;
}

/* write_colour_image:
   Write a BIMAGE object to file as a colour image
*/
void write_colour_image(
    BIMAGE * bimage_red,
    BIMAGE * bimage_green,
    BIMAGE * bimage_blue,
    const char * filename
    )
{
    /* Various objects needed by ImageMagick */
    ExceptionInfo exception;
    Image * image;
    ImageInfo * image_info;
    PixelPacket * pixel_packet;

    int i, num_pixels;

    /* Initialise ImageMagick.  This just tells it what directory (*argv) we're running from */
    MagickCoreGenesis((char *)NULL, MagickFalse);

    /* Initialise the ExceptionInfo - ImageMagick's way of implementing exceptions in C */
    GetExceptionInfo(&exception);

    /* Construct an image_info object by 'cloning' a null image_info object */
    image_info = CloneImageInfo((ImageInfo *) NULL);

    /*** Construct an ImageMagick Image of the appropriate dimensions ***/
    image = ConstituteImage(
        bimage_red->dim->buf[0],			// nx
        bimage_red->dim->buf[1],			// ny
        "I",							// Grayscale image, contains only "Intensity"
        CharPixel, 						// Data type
        bimage_red->buf,					// Image data
        &exception						// Exception object
	);

    /*** Copy pixel values from BIMAGE object to ImageMagick's Image object ***/
    num_pixels = BVECT_prod(bimage_red->dim);
    pixel_packet = GetAuthenticPixels(image, 0, 0, image->columns, image->rows, &exception);
    if (pixel_packet == NULL) {
        printf("write_grayscale_image: GetImagePixels failed!\n");
        return;
    }
    for (i = 0; i < num_pixels; ++i) {
        pixel_packet[i].red = (unsigned short int)bimage_red->buf[i];
        pixel_packet[i].green = (unsigned short int)bimage_green->buf[i];
        pixel_packet[i].blue = (unsigned short int)bimage_blue->buf[i];
    }
    /* Synchronise the buffer with the contour image!! */
    SyncAuthenticPixels(image, &exception);

    /*** Write the image to file ***/
    /* To save an image to a specific location, store the filename in the image itself */
    (void) strcpy(image->filename, filename);
    /* Actually write the image */
    // image->depth=8; // force 8-bit output
    SetImageDepth(image, 8);
    WriteImage((ImageInfo *)NULL, image);

    DestroyImage(image);

    /* Destroy the image info object */
    image_info = DestroyImageInfo(image_info);
    
    /* Finalise the ExceptionInfo */
    DestroyExceptionInfo(&exception);
    
    /* Finalise ImageMagick */
    MagickCoreTerminus();
    
    return;
}

/* normalise_contrast:
   Normalise the image contrast.  Specific to our dataset!
*/
void normalise_contrast(
    BIMAGE * image
    )
{
    int i, num_pixels;
    int min_value, max_value;
	
    num_pixels = BVECT_prod(image->dim);
	
#ifdef NORMALISE_CONSTRAST_DYNAMIC
    /* Dynamically compute range of image */
    min_value = MAX_PIXEL_VALUE;
    max_value = 0;
    for (i = 0; i < num_pixels; ++i) {
        if (image->buf[i] < min_value) min_value = image->buf[i];
        if (image->buf[i] > max_value) max_value = image->buf[i];
    }

    printf("min_value = %f\n", min_value);
    printf("max_value = %f\n", max_value);
#endif
#ifdef NORMALISE_CONTRAST_STATIC
    /* Known data range */
    min_value = MIN_MRI_PIXEL_VALUE;
    max_value = MAX_MRI_PIXEL_VALUE;
#endif

    /* Normalise range of image */
    for (i = 0; i < num_pixels; ++i) {
#ifdef NORMALISE_CONTRAST_STATIC
        if (image->buf[i] < min_value) image->buf[i] = min_value;
        if (image->buf[i] > max_value) image->buf[i] = max_value;
#endif
        image->buf[i] = MAX_PIXEL_VALUE * (image->buf[i] - min_value) / (max_value - min_value);
    }
}
