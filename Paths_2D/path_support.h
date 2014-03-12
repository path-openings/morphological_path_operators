/*
 * File:		path_support.h
 *
 * Written by:		Ben Appleton
 *
 * Date:		July 2004

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
 

/*********************************************************************************************
 path_support.h
 ------

  DESCRIPTION:
  Supporting structures and code for the grayscale path opening transform.

  HISTORY:
  Created by Ben Appleton (July 2004)
  Contact: appleton@itee.uq.edu.au
**********************************************************************************************/

#ifndef PATH_SUPPORT_H
#define PATH_SUPPORT_H

#ifndef MAX
	#define MAX(x, y) ((x) > (y) ? (x) : (y))
#endif
#ifndef MIN
	#define MIN(x, y) ((x) < (y) ? (x) : (y))
#endif

/* Specify the pixel type used by my path openings */
/* #define GPOT_PIX_TYPE int */
#define GPOT_PIX_TYPE unsigned char
#define GPOT_PIX_TYPE_NUM 256

/* pathopen now shares code with GPOT */
#define PATHOPEN_PIX_TYPE GPOT_PIX_TYPE
#define PATHOPEN_PIX_TYPE_NUM GPOT_PIX_TYPE_NUM

/* Radix sort or qsort? */
#define GPOT_RADIXSORT
/* #define GPOT_QSORT */

#define PATH_GRANULOMETRY_MIN_ALLOCATED_LENGTH 10

/*************************************** QUEUES ******************************************/
/* A unique queue for each row */
typedef struct {
	int * buffer;
	int * * index;
	int * length;
	int num_rows;
} ROW_QUEUE;

/* A working queue */
typedef struct {
	int * index;
	int length;
} QUEUE;

/******************************** PATH_GRANULOMETRY **************************************/
/*
 * Store the set of path lengths and associated thresholds for each point.  The collection
 * of these for each pixel represents the grayscale path opening transform.
 *
 */
typedef struct {
	int length;						/* Length of arrays */
	int allocated_length;			/* Allocated length of arrays */
	
	unsigned char * buf;			/* Combined memory pool for arrays */
	int * path_length;				/* Base pointer for array of path lengths */
	GPOT_PIX_TYPE * threshold;		/* Base pointer for array of associated thresholds */
} PATH_GRANULOMETRY;


/* - PATH_GRANULOMETRY_constructor:
 * Construct a granulometric curve of the specific length.  Allocates double memory.
 */
PATH_GRANULOMETRY * PATH_GRANULOMETRY_constructor(
	int length
);

/* - PATH_GRANULOMETRY_destructor:
 * Deallocate internal memory and object memory
 */
void PATH_GRANULOMETRY_destructor(
	PATH_GRANULOMETRY * this_struct
);

/* - add_point:
 * Add a point to the given granulometric curve
 */
void add_point(
	int path_length,						/* The path_length of the new point */
	GPOT_PIX_TYPE threshold,				/* The grayscale threshold of the new point */
	PATH_GRANULOMETRY * this_struct			/* The 'this' pointer */
);

/* - path_length_to_threshold:
 * Given a specified path length, determine the corresponding threshold in the granulometric curve
 */
GPOT_PIX_TYPE path_length_to_threshold(
	int path_length,						/* The desired path length */
	PATH_GRANULOMETRY * this_struct			/* The 'this' pointer */
);

/* - threshold_to_path_length:
 * Given a specified grayscale threshold, determine the corresponding path length in the granulometric curve
 */
int threshold_to_path_length(
	GPOT_PIX_TYPE threshold,				/* The desired grayscale threshold */
	PATH_GRANULOMETRY * this_struct			/* The 'this' pointer */
);

/* - PATH_GRANULOMETRY_print:
 * Print the given structure
 */
void PATH_GRANULOMETRY_print(
	PATH_GRANULOMETRY * this_struct			/* Structure to print */
);

/* - merge:
 * Merge two granulometric curves, constructing and returning a new one.
 */
PATH_GRANULOMETRY * merge(
	PATH_GRANULOMETRY * input_struct_a,			/* Structure to merge, "A" */
	PATH_GRANULOMETRY * input_struct_b			/* Structure to merge, "B" */
);


/******************************** UTILITY FUNCTION PROTOTYPES **************************************/
/* Sort an image by its pixel values (monotonic transform to 0, 1, ...) */
void image_sort(
	GPOT_PIX_TYPE * input_image,
	int num_pixels,
	int * sorted_indices
);

/* Sort an array of pointers to pixels */
int pointer_value_comparison(
	const void * a,
	const void * b
);

/* Sort an array of integer indices */
int integer_comparison(
	const void * a,
	const void * b
);

/* Transpose an image, possibly 'in-place' (internally allocates extra memory) */
void transpose_image(
	void * input_image,
	int nx, int ny,
	int num_bytes_per_element,
	void * output_image
);

/* Transform an array of pixel indices according to a transposition of the image */
void transpose_indices(
	int * input_indices,
	int nx, int ny,
	int * output_indices
);

/* Vertically flip an image, possibly 'in-place' (internally allocates extra memory) */
void flip_image(
	void * input_image,
	int nx, int ny,
	int num_bytes_per_element,
	void * output_image
);

/* Transform an array of pixel indices according to a vertical flip of the image */
void flip_indices(
	int * input_indices,
	int nx, int ny,
	int * output_indices
);

#endif // PATH_SUPPORT_H
