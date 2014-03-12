/*
 * File:		path_support.c
 *

 Copyright Benjamin Appleton and Hugues Talbot, July 2004

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
 gpot.c
 ------

  DESCRIPTION:
  A grayscale path opening transform.

  REFERENCE:
  Previous work: Hugues Talbot, Michael Buckley, Henk Heijmans (in no particular order)

  HISTORY:
  Modified from pathopen.c July 2004.  Both by Ben Appleton.
  Contact: appleton@itee.uq.edu.au
**********************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "path_support.h"

/**************************** PATH_GRANULOMETRY IMPLEMENTATIONS *******************************/
/* - PATH_GRANULOMETRY_constructor:
 * Construct a granulometric curve of the specific length.  Allocates double memory.
 */
PATH_GRANULOMETRY * PATH_GRANULOMETRY_constructor(
	int length
)
{
	PATH_GRANULOMETRY * this_struct;

	/* Allocate memory for the structure itself */
	this_struct = (PATH_GRANULOMETRY *)malloc(sizeof(PATH_GRANULOMETRY));

	/* Allocate the internal memory */
	this_struct->length = length >= 0 ? length : 0;
	this_struct->allocated_length = MAX(this_struct->length * 2, PATH_GRANULOMETRY_MIN_ALLOCATED_LENGTH);

	/* Keep both arrays (path_length, threshold) inside single larger array */
	this_struct->buf = (void *)malloc(this_struct->allocated_length * (sizeof(int) + sizeof(GPOT_PIX_TYPE)));
	this_struct->path_length = (int *)this_struct->buf;
	this_struct->threshold = (GPOT_PIX_TYPE *)(this_struct->buf + this_struct->allocated_length * sizeof(int));

	return this_struct;
}

/* - PATH_GRANULOMETRY_destructor:
 * Deallocate internal memory and object memory
 */
void PATH_GRANULOMETRY_destructor(
	PATH_GRANULOMETRY * this_struct
)
{
	free((void *)this_struct->buf);
	free((void *)this_struct);
}

/* - add_point:
 * Add a point to the given granulometric curve, eliding into existing points if appropriate
 */
void add_point(
	int path_length,						/* The path_length of the new point */
	GPOT_PIX_TYPE threshold,				/* The grayscale threshold of the new point */
	PATH_GRANULOMETRY * this_struct			/* The 'this' pointer */
)
{
	/* As points are passed in, lengths monotonically non-increase, thresholds monotonically non-decrease */
	/* Check if a new point needs to be added to the granulometric curve */
	if (this_struct->length > 0) {
		/* == length, > threshold -> update existing point */
		if (path_length == this_struct->path_length[this_struct->length - 1]) {
			if (threshold > this_struct->threshold[this_struct->length - 1]) {
				this_struct->threshold[this_struct->length - 1] = threshold;
			}
			return;
		}
		/* == threshold -> skip */
		if (threshold == this_struct->threshold[this_struct->length - 1]) {
			return;
		}
	}

	/* Create a new point in the granulometric curve */
	this_struct->path_length[this_struct->length] = path_length;
	this_struct->threshold[this_struct->length] = threshold;
	++(this_struct->length);

	/* Check if we need to allocate more memory */
	if (this_struct->length == this_struct->allocated_length) {
		PATH_GRANULOMETRY * new_struct;

		/* Construct a new structure with the same length but twice the allocated memory */
		new_struct = PATH_GRANULOMETRY_constructor(this_struct->length);

		/* Copy the data across from old to new */
		memcpy(new_struct->path_length, this_struct->path_length, this_struct->length * sizeof(int));
		memcpy(new_struct->threshold, this_struct->threshold, this_struct->length * sizeof(GPOT_PIX_TYPE));
		this_struct->length = new_struct->length;
		this_struct->allocated_length = new_struct->allocated_length;

		/* Now swap internal arrays, gutting the new_struct object, and delete the old structure */
		{
			void * swap = this_struct->buf;
			this_struct->buf = new_struct->buf;
			free((void *)swap);
		}
		this_struct->path_length = new_struct->path_length;
		this_struct->threshold = new_struct->threshold;
		
		/* Clean up by deleting the skeleton of the new structure */
		free((void *)new_struct);
	}
}

/* - path_length_to_threshold:
 * Given a specified path length, determine the corresponding threshold in the granulometric curve
 */
GPOT_PIX_TYPE path_length_to_threshold(
	int path_length,						/* The desired path length */
	PATH_GRANULOMETRY * this_struct			/* The 'this' pointer */
)
{
	int i;
	GPOT_PIX_TYPE value;

	if (this_struct->length == 0) {
		return (GPOT_PIX_TYPE)0;
	}

	/* Granulometric curves are stored in order of increasing threshold, decreasing path length */
	value = 0;
	for (i = 0; i < this_struct->length; ++i) {
		/* Once we find the correct segment of the granulometric curve, halt */
		if (this_struct->path_length[i] < path_length) break;

		value = this_struct->threshold[i];
	}

	return value;
}

/* - threshold_to_path_length:
 * Given a specified grayscale threshold, determine the corresponding path length in the granulometric curve
 */
int threshold_to_path_length(
	GPOT_PIX_TYPE threshold,				/* The desired grayscale threshold */
	PATH_GRANULOMETRY * this_struct			/* The 'this' pointer */
)
{
	int i;
	int path_length;

	if (this_struct->length == 0) {
		return (int)0;
	}

	/* Granulometric curves are stored in order of increasing threshold */
	path_length = 0;
	for (i = 0; i < this_struct->length; ++i) {
		path_length = this_struct->path_length[i];

		/* Once we find the correct segment of the granulometric curve, halt */
		if (this_struct->threshold[i] > threshold) break;
	}

	return path_length;
}

/* - PATH_GRANULOMETRY_print:
 * Print the given structure
 */
void PATH_GRANULOMETRY_print(
	PATH_GRANULOMETRY * this_struct			/* Structure to print */
)
{
	int i;

	printf("[");
	for (i = 0; i < this_struct->length; ++i) {
		printf("(%i, %i)", this_struct->path_length[i], (int)this_struct->threshold[i]);
	}
	printf("]");
}

/* - merge:
 * Merge two granulometric curves, constructing and returning a new one.
 */
PATH_GRANULOMETRY * merge(
	PATH_GRANULOMETRY * input_struct_a,			/* Structure to merge, "A" */
	PATH_GRANULOMETRY * input_struct_b			/* Structure to merge, "B" */
)
{
	int i;
	int index_a, index_b, index_output;
	int cur_path_length;
	int num_points;
	int cur_threshold;		/* This is slightly naughty (GPOT_PIX_TYPE could be float), but gets around more serious problems */
	PATH_GRANULOMETRY * this_struct;

	/* Allocate structure which will be returned.  Initially empty */
	this_struct = PATH_GRANULOMETRY_constructor(0);

	/* A conceptual algorithm to merge the two lists:
	 * Put all points from both lists into the one list
	 * If any point has simultaneously <= path_length and <= threshold than some other point,
	 * we can remove that point.  All other points must remain.
	 */

	/* Note - Points are in increasing order of threshold due to the original code for pathopen.
	 * So, they are also in decreasing order of length.
	 */
	index_a = index_b = index_output = 0;
	cur_path_length = (1 << 30);
	cur_threshold = -1;

	/* Set up dummy endpoints to reduce if statements in inner loop */
	/* Note - PATH_GRANULOMETRY structure always has a blank node allocated at end of list, use that */
	input_struct_a->path_length[input_struct_a->length] = -1;
	input_struct_b->path_length[input_struct_b->length] = -1;
	num_points = input_struct_a->length + input_struct_b->length;
	for (i = 0; i < num_points; ++i) {
		int new_path_length;
		GPOT_PIX_TYPE new_threshold;

		/* Competitively consume the lists in order of length (descending order) */
		if (input_struct_a->path_length[index_a] > input_struct_b->path_length[index_b]) {
			new_path_length = input_struct_a->path_length[index_a];
			new_threshold = input_struct_a->threshold[index_a];
			++index_a;
		} else {
			new_path_length = input_struct_b->path_length[index_b];
			new_threshold = input_struct_b->threshold[index_b];
			++index_b;
		}

		/* Can this point add a new output point? */
		if (new_threshold > cur_threshold) {
			/* If new point has lower length and greater threshold, add point */
			if (new_path_length < cur_path_length) {
				/* Add a new output point */
				add_point(new_path_length, new_threshold, this_struct);
				cur_path_length = new_path_length;
				cur_threshold = (int)new_threshold;
			} else {
				/* Update the threshold of the existing output point */
				this_struct->threshold[this_struct->length - 1] = new_threshold;
				cur_threshold = (int)new_threshold;
			}
		}
	}

	return this_struct;
}


/****************************** UTILITY FUNCTION IMPLEMENTATIONS ***********************************/
/* Compare two pointers to image pixels (sorts ascending) */
int pointer_value_comparison(
	const void * a,
	const void * b
)
{
	return (
		(int)(**((GPOT_PIX_TYPE **)a)) - (int)(**((GPOT_PIX_TYPE **)b))
	);
}


/* Sort an array of ints (ascending). */
int integer_comparison(
	const void * a,
	const void * b
)
{
	return (
		*((int *)a) - *((int *)b)
	);
}


/* Sort an image by its pixel values (monotonic transform to 0, 1, ...) */
void image_sort(
	GPOT_PIX_TYPE * input_image,
	int num_pixels,
	int * sorted_indices
)
{
	int i;

#ifdef GPOT_QSORT
	GPOT_PIX_TYPE * * sorted_pointers;

	/* Allocate memory */
	sorted_pointers = (GPOT_PIX_TYPE **)malloc(num_pixels * sizeof(GPOT_PIX_TYPE *));

	/* Create an array of pointers to sort */
	for (i = 0; i < num_pixels; ++i) {
		sorted_pointers[i] = input_image + i;
	}

	/* Sort the array of pointers.  WARNING: qsort's implementation is system dependent, may be slow! */
	qsort(sorted_pointers, num_pixels, sizeof(GPOT_PIX_TYPE *), pointer_value_comparison);

	/* Convert pointers to indices */
	for (i = 0; i < num_pixels; ++i) {
		sorted_indices[i] = (int)(sorted_pointers[i] - input_image);
	}

	/* Free memory */
	free((void *)sorted_pointers);
#endif
#ifdef GPOT_RADIXSORT
	/* Store the length of each run of pixel values */
	int length[GPOT_PIX_TYPE_NUM];
	int * data[GPOT_PIX_TYPE_NUM];

	/* How many pixels of each value? */
	/* Note: Here I count, allocate, and count again
	 - slightly wasteful, but not enough to bother recoding */
	memset(length, 0, GPOT_PIX_TYPE_NUM * sizeof(int));
	for (i = 0; i < num_pixels; ++i) {
		length[input_image[i]]++;
	}

	/* Now set up data pointers into sorted_indices block */
	data[0] = sorted_indices;
	for (i = 1; i < GPOT_PIX_TYPE_NUM; ++i) {
		data[i] = data[i - 1] + length[i - 1];
	}

	/* Now sort the image, destroying the data pointers which are not used after this point */
	for (i = 0; i < num_pixels; ++i) {
		*(data[input_image[i]]) = i;
		++data[input_image[i]];
	}
#endif
}


/*
	Transpose an image, possibly 'in place'
*/
void transpose_image(
	void * input_image,
	int nx, int ny,
	int num_bytes_per_element,
	void * output_image
)
{
	char in_place;
	void * temporary_image;
	int i, num_pixels;

	num_pixels = nx * ny;

	in_place = (input_image == output_image);

	/* Allocate memory */
	if (in_place) {
		temporary_image = (void *)malloc(num_pixels * num_bytes_per_element);
	} else {
		temporary_image = output_image;
	}

	/* Transpose from input image to temporary image */
	/* Templates would be nifty here! */
#ifdef GPOT_USE_MEMCPY
	for (i = 0; i < num_pixels; ++i) {
		int x, y, new_index;

		/* Compute destination address */
		x = i % nx;
		y = i / nx;
		new_index = y + ny * x;

		memcpy(temporary_image + new_index * num_bytes_per_element, input_image + i * num_bytes_per_element, num_bytes_per_element);
	}
#else
	switch(num_bytes_per_element) {
		case 1:
			for (i = 0; i < num_pixels; ++i) {
				int x, y, new_index;

				/* Compute destination address */
				x = i % nx;
				y = i / nx;
				new_index = y + ny * x;

				((char *)temporary_image)[new_index] = ((char *)input_image)[i];
			}
			break;
		case 2:
			for (i = 0; i < num_pixels; ++i) {
				int x, y, new_index;

				/* Compute destination address */
				x = i % nx;
				y = i / nx;
				new_index = y + ny * x;

				((short int *)temporary_image)[new_index] = ((short int *)input_image)[i];
			}
			break;
		case 4:
			for (i = 0; i < num_pixels; ++i) {
				int x, y, new_index;

				/* Compute destination address */
				x = i % nx;
				y = i / nx;
				new_index = y + ny * x;

				((long int *)temporary_image)[new_index] = ((long int *)input_image)[i];
			}
			break;
		default:
			for (i = 0; i < num_pixels; ++i) {
				int x, y, new_index;

				/* Compute destination address */
				x = i % nx;
				y = i / nx;
				new_index = y + ny * x;

				memcpy(
					(unsigned char *)(temporary_image) + new_index * num_bytes_per_element, 
					(unsigned char *)(input_image) + i * num_bytes_per_element, 
					num_bytes_per_element
				);
			}
			break;
	}
#endif

	/* Free memory */
	if (in_place) {
		memcpy(output_image, temporary_image, num_pixels * num_bytes_per_element);
		free((void *)temporary_image);
	}
}


/* - transpose_indices:
	Transform an array of pixel indices according to a transposition.  Doesn't
	care about in-place/out-of-place
*/
void transpose_indices(
	int * input_indices,
	int nx, int ny,
	int * output_indices
)
{
	int i, num_pixels;

	num_pixels = nx * ny;

	/* Transform each index individually */
	for (i = 0; i < num_pixels; ++i) {
		int x, y, old_index, new_index;

		/* Lookup old image index */
		old_index = input_indices[i];

		/* Compute destination address */
		x = old_index % nx;
		y = old_index / nx;
		new_index = y + ny * x;

		/* Store transposed image index */
		output_indices[i] = new_index;
	}
}


/*
	Flip an image, possibly 'in place'
*/
void flip_image(
	void * input_image,
	int nx, int ny,
	int num_bytes_per_element,
	void * output_image
)
{
	char in_place;
	void * temporary_image;
	int y, num_pixels;

	num_pixels = nx * ny;

	in_place = (input_image == output_image);

	/* Allocate memory */
	if (in_place) {
		temporary_image = (void *)malloc(num_pixels * num_bytes_per_element);
	} else {
		temporary_image = output_image;
	}

	/* Flip the image */
	for (y = 0; y < ny; ++y) {
		int row_base_index = nx * y;
		int new_row_base_index = nx * (ny - 1 - y);

		memcpy(
			(unsigned char *)temporary_image + new_row_base_index * num_bytes_per_element,
			(unsigned char *)input_image + row_base_index * num_bytes_per_element,
			nx * num_bytes_per_element
		);
	}

	/* Free memory */
	if (in_place) {
		memcpy(output_image, temporary_image, num_pixels * num_bytes_per_element);
		free((void *)temporary_image);
	}
}


/* - flip_indices:
	Transform an array of pixel indices according to a vertical flip.  Doesn't
	care about in-place/out-of-place
*/
void flip_indices(
	int * input_indices,
	int nx, int ny,
	int * output_indices
)
{
	int i, num_pixels;

	num_pixels = nx * ny;

	/* Transform each index individually */
	for (i = 0; i < num_pixels; ++i) {
		int x, y, old_index, new_index;

		/* Lookup old image index */
		old_index = input_indices[i];

		/* Compute destination address */
		x = old_index % nx;
		y = old_index / nx;
		new_index = x + nx * (ny - 1 - y);

		/* Store transposed image index */
		output_indices[i] = new_index;
	}
}


