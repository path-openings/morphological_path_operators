/*
 * File:		pathopen.cxx
 *
  Copyright Benjamin Appleton and Hugues Talbot, April 2004

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
 pathopen.cxx
 ------

  DESCRIPTION:
  Path opening for grayscale images.

  REFERENCE:
  Previous work: Hugues Talbot, Michael Buckley, Henk Heijmans (in no particular order)

  HISTORY:
  Modified by Ben Appleton (July 2005):
	- Performs incomplete path openings 
	- Converted trivially from C to C++
  Created by Ben Appleton (February 2004)
  Contact: appleton@itee.uq.edu.au
**********************************************************************************************/

#pragma warning(disable:4786)
//#pragma warning(disable: 4237 4284 4290 4514 4786) // Disable silly 'warning's about things like long templated variable names

//#define DEBUGGING

#include <conio.h>

#include "pathopen.h"
#include "path_queue.h"

#include <algorithm>
#include <iostream>
using namespace std;

/* - pathopen:
	Perform a path opening on an image.  Main interface, calls all subfunctions.
*/
int pathopen(
	PATHOPEN_PIX_TYPE * input_image,				/* The input image */
	int nx, int ny,									/* Image dimensions */
	int L,											/* The threshold line length */
	int K,											/* The maximum number of gaps in the path */
	PATHOPEN_PIX_TYPE * output_image				/* Output image */
)
{
	int i, num_pixels;

	PATHOPEN_PIX_TYPE * accumulator_image;

	PATHOPEN_PIX_TYPE * transposed_input_image;
	PATHOPEN_PIX_TYPE * flipped_input_image;

	int * sorted_indices;
	int * transposed_sorted_indices;
	int * flipped_sorted_indices;

	num_pixels = nx * ny;

	/* Allocate memory */
	accumulator_image = (PATHOPEN_PIX_TYPE *)malloc(num_pixels * sizeof(PATHOPEN_PIX_TYPE));

	transposed_input_image = (PATHOPEN_PIX_TYPE *)malloc(num_pixels * sizeof(PATHOPEN_PIX_TYPE));
	flipped_input_image = (PATHOPEN_PIX_TYPE *)malloc(num_pixels * sizeof(PATHOPEN_PIX_TYPE));

	sorted_indices = (int *)malloc(num_pixels * sizeof(int));
	transposed_sorted_indices = (int *)malloc(num_pixels * sizeof(int));
	flipped_sorted_indices = (int *)malloc(num_pixels * sizeof(int));

	/* Sort the image pixels, storing pixel indices */
	image_sort(input_image, num_pixels, sorted_indices);

	/* Create a transposed copy of the original image */
	transpose_image((void *)input_image, nx, ny, sizeof(PATHOPEN_PIX_TYPE), (void *)transposed_input_image);
	transpose_indices(sorted_indices, nx, ny, transposed_sorted_indices);

	/* Create a flipped copy of the original image */
	flip_image((void *)input_image, nx, ny, sizeof(PATHOPEN_PIX_TYPE), (void *)flipped_input_image);
	flip_indices(sorted_indices, nx, ny, flipped_sorted_indices);

	/* Vertical path opening */
	vert_pathopen(input_image, sorted_indices, nx, ny, L, K, output_image);

	/* ++diagonal path opening */
	diag_pathopen(input_image, sorted_indices, nx, ny, L, K, accumulator_image);

	/* Accumulate results */
	for (i = 0; i < num_pixels; ++i) {
		output_image[i] = MAX(output_image[i], accumulator_image[i]);
	}

	/* Horizontal path opening */
	vert_pathopen(transposed_input_image, transposed_sorted_indices, ny, nx, L, K, accumulator_image);
	transpose_image((void *)accumulator_image, ny, nx, sizeof(PATHOPEN_PIX_TYPE), (void *)accumulator_image);
	/* Accumulate into output */
	for (i = 0; i < num_pixels; ++i) {
		output_image[i] = MAX(output_image[i], accumulator_image[i]);
	}

	/* +-diagonal path opening */
	diag_pathopen(flipped_input_image, flipped_sorted_indices, nx, ny, L, K, accumulator_image);
	flip_image((void *)accumulator_image, nx, ny, sizeof(PATHOPEN_PIX_TYPE), (void *)accumulator_image);
	/* Accumulate into output */
	for (i = 0; i < num_pixels; ++i) {
		output_image[i] = MAX(output_image[i], accumulator_image[i]);
	}

	/* Free allocated memory */
	free((void *)sorted_indices);
	free((void *)transposed_sorted_indices);
	free((void *)flipped_sorted_indices);
	free((void *)accumulator_image);
	free((void *)transposed_input_image);
	free((void *)flipped_input_image);

	return 0;
}


/* A path opening in the vertical direction.
	Conjugate with transpose to perform horizontal path openings
*/
static int vert_pathopen(
	PATHOPEN_PIX_TYPE * input_image,					/* The input image */
	int * sorted_indices,								/* Monotonic transform to [0, 1, ...] of input image */
	int nx, int ny,										/* Image dimensions */
	int L,												/* The threshold line length */
	int K,												/* The maximum gap number */
	PATHOPEN_PIX_TYPE * output_image					/* Output image */
)
{
	int i, k, x, y, index, new_index, sort_index, num_pixels;

	/************************************** Allocation **********************************************/
	num_pixels = nx * ny;
	int nk = K + 1;

	/* Construct queueing system */
	Path_Queue path_queue_up(nk, ny, nx);
	Path_Queue path_queue_down(nk, ny, nx);

	/* Dynamic binary input image */
	char * bin_input_image = (char *)malloc(num_pixels * sizeof(char));

	/* in_queue flags */
	char * in_queue_up = (char *)malloc(num_pixels * nk * sizeof(char));
	char * in_queue_down = (char *)malloc(num_pixels * nk * sizeof(char));

	/* Chain length images [k + nk * pixel_index].  These don't include the current pixel. */
	int * chain_image_up = (int *)malloc(num_pixels * nk * sizeof(int));
	int * chain_image_down = (int *)malloc(num_pixels * nk * sizeof(int));

	// At each pixel, we store the vector of binary outputs indexed by gap number of upward chain
	char * bin_output_image_array = (char *)malloc(num_pixels * nk * sizeof(char));
	// Also count the vector of binary outputs, to note when they are all extinguished (boolean PQ!)
	char * bin_output_image_count = (char *)malloc(num_pixels * sizeof(char));

	/************************************** Initialisation **********************************************/
	/* Dynamic binary threshold image is initially all 1's */
	memset(bin_input_image, 1, num_pixels * sizeof(char));

	/* Queue initially empty */
	memset(in_queue_up, 0, num_pixels * nk * sizeof(char));
	memset(in_queue_down, 0, num_pixels * nk * sizeof(char));

	/* Initialise the chain lengths */
	for (y = 0; y < ny; ++y) {
		int up_length = y;
		int down_length = (ny - 1) - y;

#ifdef PATHOPEN_LENGTH_HEURISTIC
		/* Threshold at L - 1 */
		if (up_length > L - 1) up_length = L - 1;
		if (down_length > L - 1) down_length = L - 1;
#endif

		/* Memset */
		for (x = 0, index = nx * y; x < nx; ++x, ++index) 
			for (k = 0; k < nk; ++k)
				chain_image_up[k + nk * index] = up_length;
		for (x = 0, index = nx * y; x < nx; ++x, ++index) 
			for (k = 0; k < nk; ++k)
				chain_image_down[k + nk * index] = down_length;
	}

	/* Binary output vector at each pixel */
	memset(bin_output_image_array, 1, num_pixels * nk * sizeof(char));
	memset(bin_output_image_count, nk, num_pixels * sizeof(char));

	/* Set output image to default value (0 here) */
	memset(output_image, 0, num_pixels * sizeof(PATHOPEN_PIX_TYPE));
	/****************************************************************************************************/

	/* For each threshold from smallest to largest */
	sort_index = 0;
	while(sort_index < num_pixels) {
		PATHOPEN_PIX_TYPE threshold;

		/*********************************** Process threshold pixels *****************************************/
		threshold = input_image[sorted_indices[sort_index]];
#ifdef DEBUGGING
		cout << "Threshold = " << (int)threshold << endl;
#endif
		while(input_image[sorted_indices[sort_index]] == threshold) {
			/* Collect into rows for enqueueing */
			vector< vector<PIXEL_INDEX_TYPE> > new_row_queue_down(nk);
			vector< vector<PIXEL_INDEX_TYPE> > new_row_queue_up(nk);
			int row_y = sorted_indices[sort_index] / nx;
#ifdef DEBUGGING
			cout << "y = " << row_y << endl;
#endif

			while(input_image[sorted_indices[sort_index]] == threshold && sorted_indices[sort_index] / nx == row_y) {
				/* Extract index and coordinates */
				index = sorted_indices[sort_index];
				y = index / nx;
				x = index % nx;

#ifdef DEBUGGING
				cout << "\tx = " << x << endl;
#endif

				/* Directly perform the changes to this threshold pixel */
				if (bin_input_image[index]) {
#ifdef DEBUGGING
					cout << "\tRemoving..." << endl;
					cout << "\tcount = " << (int)bin_output_image_count[index] << endl;
#endif
					// Remove this pixel from the binary input image
					bin_input_image[index] = 0;

					// Update the outputs
					if (bin_output_image_count[index] > 0) {
						// Update the output flags for each gap index (and count)
						bin_output_image_count[index] = 0;
						for (k = 0; k < K; ++k) {
							bin_output_image_array[k + nk * index] = 
								(chain_image_up[k + nk * index] + chain_image_down[(K - 1 - k) + nk * index] + 1) >= L;
							bin_output_image_count[index] += bin_output_image_array[k + nk * index];
						}

						// If all paths have been extinguished, update output
						if (bin_output_image_count[index] == 0) {
							output_image[index] = threshold;
						}
					}

					/* Enqueue downward pixels for update */
					if (y < ny - 1) {
						if (x > 0) {
							new_index = index + nx - 1;
							
							// Enqueue this pixel for all k
							for (k = 0; k < nk; ++k) {
								if (!in_queue_down[k + nk * new_index]) {
									in_queue_down[k + nk * new_index] = 1;
									new_row_queue_down[k].push_back(x - 1);
								}
							}
						}
						new_index = index + nx;
						// Enqueue this pixel for all k
						for (k = 0; k < nk; ++k) {
							if (!in_queue_down[k + nk * new_index]) {
								in_queue_down[k + nk * new_index] = 1;
								new_row_queue_down[k].push_back(x);
							}
						}
						if (x < nx - 1) {
							new_index = index + nx + 1;
							
							// Enqueue this pixel for all k
							for (k = 0; k < nk; ++k) {
								if (!in_queue_down[k + nk * new_index]) {
									in_queue_down[k + nk * new_index] = 1;
									new_row_queue_down[k].push_back(x + 1);
								}
							}
						}
					}

					/* Enqueue upward pixels for update */
					if (y > 0) {
						if (x > 0) {
							new_index = index - nx - 1;
							
							// Enqueue this pixel for all k
							for (k = 0; k < nk; ++k) {
								if (!in_queue_up[k + nk * new_index]) {
									in_queue_up[k + nk * new_index] = 1;
									new_row_queue_up[k].push_back(x - 1);
								}
							}
						}
						new_index = index - nx;
						// Enqueue this pixel for all k
						for (k = 0; k < nk; ++k) {
							if (!in_queue_up[k + nk * new_index]) {
								in_queue_up[k + nk * new_index] = 1;
								new_row_queue_up[k].push_back(x);
							}
						}
						if (x < nx - 1) {
							new_index = index - nx + 1;
							
							// Enqueue this pixel for all k
							for (k = 0; k < nk; ++k) {
								if (!in_queue_up[k + nk * new_index]) {
									in_queue_up[k + nk * new_index] = 1;
									new_row_queue_up[k].push_back(x + 1);
								}
							}
						}
					}
				}
				
				// Next pixel in row			
				++sort_index;
				if (sort_index >= num_pixels) break;
			}
#ifdef DEBUGGING
			cout << "Merging row queues" << endl;
#endif

			// Merge row queues into path queues
			if (row_y + 1 < ny) {
#ifdef DEBUGGING
				cout << "Down queues" << endl;
#endif
				for (k = 0; k < nk; ++k) {
					if (new_row_queue_down[k].size() > 0) {
#ifdef DEBUGGING
						cout << "k = " << k << ":";
						for (i = 0; i < new_row_queue_down[k].size(); ++i) {
							cout << new_row_queue_down[k][i] << " ";
						}
						cout << endl;
#endif
						path_queue_down.merge_row(new_row_queue_down[k], k, row_y + 1);
					}
				}
			}

			if (row_y - 1 >= 0) {
#ifdef DEBUGGING
				cout << "Up queues" << endl;
#endif
				for (k = 0; k < nk; ++k) {
					if (new_row_queue_up[k].size() > 0) {
#ifdef DEBUGGING
						cout << "k = " << k << ":";
						for (i = 0; i < new_row_queue_up[k].size(); ++i) {
							cout << new_row_queue_up[k][i] << " ";
						}
						cout << endl;
#endif
						path_queue_up.merge_row(new_row_queue_up[k], k, row_y - 1);
					}
				}
			}

#ifdef DEBUGGING
			cout << endl;
#endif
			if (sort_index >= num_pixels) break;
		}

		/*************************************** Downward sweep *********************************************/
		/* Propagate changes at current threshold down the image */
#ifdef DEBUGGING
		cout << "DOWNWARD SWEEP - before" << endl;
		dump_state(
			nx,									// Image dimensions etc.
			ny,
			nk,
			path_queue_up,						// Queueing structures
			path_queue_down,
			in_queue_up,
			in_queue_down,
			bin_input_image,					// Input/output binary images
			bin_output_image_array,
			bin_output_image_count,
			chain_image_up,						// Up/down chain lengths
			chain_image_down
		);
#endif

		for (k = 0; k < nk; ++k) {
#ifdef DEBUGGING
			cout << "k = " << k << endl;
#endif
			for (y = 1; y < ny; ++y) {
				vector<PIXEL_INDEX_TYPE> & row_queue = path_queue_down.q[k][y];
				if (row_queue.size() == 0) continue;
#ifdef DEBUGGING
				cout << "\ty = " << y << endl;
#endif

				vector<PIXEL_INDEX_TYPE> new_row_queue_cur_k;
				vector<PIXEL_INDEX_TYPE> new_row_queue_next_k;

				/* Perform updates on points in row_queue, propagating changes to the new row queues */
				for (i = 0; i < row_queue.size(); ++i) {
					/* Extract x-coordinate and pixel index */
					x = row_queue[i];
					index = x + nx * y;

#ifdef DEBUGGING
					cout << "\t\tx = " << x << endl;
#endif

					/* Unflag -> no longer in queue */
					in_queue_down[k + nk * index] = 0;

					/* Update chain length from upward neighbours */
					// Note: Only y > 0 may be 'updated', so we are assured of the existence of previous neighbours!
					int max_prev = -1;
					// Previous level - accept a gap
					if (k > 0) {
						if (x > 0) {
							new_index = index - nx - 1;
							if (chain_image_up[k - 1 + nk * new_index] > max_prev) {
								max_prev = chain_image_up[k - 1 + nk * new_index];
							}
						}
						new_index = index - nx;
						if (chain_image_up[k - 1 + nk * new_index] > max_prev) {
							max_prev = chain_image_up[k - 1 + nk * new_index];
						}
						if (x < nx - 1) {
							new_index = index - nx + 1;
							if (chain_image_up[k - 1 + nk * new_index] > max_prev) {
								max_prev = chain_image_up[k - 1 + nk * new_index];
							}
						}
					}
					// Current level - no gap allowed
					if (x > 0) {
						new_index = index - nx - 1;
						if (bin_input_image[new_index] == 1 && chain_image_up[k + nk * new_index] > max_prev) {
							max_prev = chain_image_up[k + nk * new_index];
						}
					}
					new_index = index - nx;
					if (bin_input_image[new_index] == 1 && chain_image_up[k + nk * new_index] > max_prev) {
						max_prev = chain_image_up[k + nk * new_index];
					}
					if (x < nx - 1) {
						new_index = index - nx + 1;
						if (bin_input_image[new_index] == 1 && chain_image_up[k + nk * new_index] > max_prev) {
							max_prev = chain_image_up[k + nk * new_index];
						}
					}

					/* Update chain length? */
					if (max_prev + 1 < chain_image_up[k + nk * index]) {
#ifdef DEBUGGING
						cout << "New chain length is " << chain_image_up[k + nk * index] << endl;
#endif
						// Update chain length
						chain_image_up[k + nk * index] = max_prev + 1;

						// Propagate changes to output
						if (bin_input_image[index]) {
							char new_bin_output_flag = 
								(chain_image_up[k + nk * index] + chain_image_down[(K - k) + nk * index] + 1 >= L);
							// Did we cross the threshold?
							if (bin_output_image_array[k + nk * index] && !new_bin_output_flag) {
								// Clear the flag
								bin_output_image_array[k + nk * index] = 0;
								--bin_output_image_count[index];
								// Did this extinguish the last path?
								if (bin_output_image_count[index] == 0) {
									// Write to output
									output_image[index] = threshold;
								}
							}
						} else {
							if (K - 1 - k >= 0) {
								char new_bin_output_flag = 
									(chain_image_up[k + nk * index] + chain_image_down[(K - 1 - k) + nk * index] + 1 >= L);
								// Did we cross the threshold?
								if (bin_output_image_array[k + nk * index] && !new_bin_output_flag) {
									// Clear the flag
									bin_output_image_array[k + nk * index] = 0;
									--bin_output_image_count[index];
									// Did this extinguish the last path?
									if (bin_output_image_count[index] == 0) {
										// Write to output
										output_image[index] = threshold;
									}
								}
							}
							// Else, this pixel has already been removed
						}

						/* Propagate changes by enqueueing downward neighbours */
						if (y < ny - 1) {
							// Same layer
							if (x > 0) {
								new_index = index + nx - 1;
								if (!in_queue_down[k + nk * new_index]) {
									new_row_queue_cur_k.push_back(x - 1);
									in_queue_down[k + nk * new_index] = 1;
								}
							}
							new_index = index + nx;
							if (!in_queue_down[k + nk * new_index]) {
								new_row_queue_cur_k.push_back(x);
								in_queue_down[k + nk * new_index] = 1;
							}
							if (x < nx - 1) {
								new_index = index + nx + 1;
								if (!in_queue_down[k + nk * new_index]) {
									new_row_queue_cur_k.push_back(x + 1);
									in_queue_down[k + nk * new_index] = 1;
								}
							}
							// Down one layer
							if (k < K) {
								if (x > 0) {
									new_index = index + nx - 1;
									if (!in_queue_down[k + 1 + nk * new_index]) {
										new_row_queue_next_k.push_back(x - 1);
										in_queue_down[k + 1 + nk * new_index] = 1;
									}
								}
								new_index = index + nx;
								if (!in_queue_down[k + 1 + nk * new_index]) {
									new_row_queue_next_k.push_back(x);
									in_queue_down[k + 1 + nk * new_index] = 1;
								}
								if (x < nx - 1) {
									new_index = index + nx + 1;
									if (!in_queue_down[k + 1 + nk * new_index]) {
										new_row_queue_next_k.push_back(x + 1);
										in_queue_down[k + 1 + nk * new_index] = 1;
									}
								}
							}
						}
					}
				}
				// Wipe old queue
				row_queue.resize(0);

				/* Merge new row queues into existing queues */
				if (y + 1 < ny) {
					if (new_row_queue_cur_k.size() > 0) {
						path_queue_down.merge_row(new_row_queue_cur_k, k, y + 1);
					}
					if (new_row_queue_next_k.size() > 0) {
						path_queue_down.merge_row(new_row_queue_next_k, k + 1, y + 1);
					}
				}
			}
		}
#ifdef DEBUGGING
		cout << "DOWNWARD SWEEP - after" << endl;
		dump_state(
			nx,									// Image dimensions etc.
			ny,
			nk,
			path_queue_up,						// Queueing structures
			path_queue_down,
			in_queue_up,
			in_queue_down,
			bin_input_image,					// Input/output binary images
			bin_output_image_array,
			bin_output_image_count,
			chain_image_up,						// Up/down chain lengths
			chain_image_down
		);
#endif

		/*************************************** Upward sweep *********************************************/
#ifdef DEBUGGING
		cout << "UPWARD SWEEP - before" << endl;
		dump_state(
			nx,									// Image dimensions etc.
			ny,
			nk,
			path_queue_up,						// Queueing structures
			path_queue_down,
			in_queue_up,
			in_queue_down,
			bin_input_image,					// Input/output binary images
			bin_output_image_array,
			bin_output_image_count,
			chain_image_up,						// Up/down chain lengths
			chain_image_down
		);
#endif

		/* Propagate changes at current threshold up the image */
		for (k = 0; k < nk; ++k) {
#ifdef DEBUGGING
			cout << "k = " << k << endl;
#endif
			for (y = ny - 2; y >= 0; --y) {
				vector<PIXEL_INDEX_TYPE> & row_queue = path_queue_up.q[k][y];
				if (row_queue.size() == 0) continue;

#ifdef DEBUGGING
				cout << "\ty = " << y << endl;
#endif

				vector<PIXEL_INDEX_TYPE> new_row_queue_cur_k;
				vector<PIXEL_INDEX_TYPE> new_row_queue_next_k;

				/* Perform updates on points in row_queue, propagating changes to the new row queues */
				for (i = 0; i < row_queue.size(); ++i) {
					/* Extract x-coordinate and pixel index */
					x = row_queue[i];
					index = x + nx * y;

#ifdef DEBUGGING
					cout << "\t\tx = " << x << endl;
#endif

					/* Unflag -> no longer in queue */
					in_queue_up[k + nk * index] = 0;

					/* Update chain length from downward neighbours */
					// Note: Only y < ny - 1 may be 'updated', so we are assured of the existence of previous neighbours!
					int max_prev = -1;
					// Previous level - accept a gap
					if (k > 0) {
						if (x > 0) {
							new_index = index + nx - 1;
							if (chain_image_down[k - 1 + nk * new_index] > max_prev) {
								max_prev = chain_image_down[k - 1 + nk * new_index];
							}
						}
						new_index = index + nx;
						if (chain_image_down[k - 1 + nk * new_index] > max_prev) {
							max_prev = chain_image_down[k - 1 + nk * new_index];
						}
						if (x < nx - 1) {
							new_index = index + nx + 1;
							if (chain_image_down[k - 1 + nk * new_index] > max_prev) {
								max_prev = chain_image_down[k - 1 + nk * new_index];
							}
						}
					}
					// Current level - no gap allowed
					if (x > 0) {
						new_index = index + nx - 1;
						if (bin_input_image[new_index] == 1 && chain_image_down[k + nk * new_index] > max_prev) {
							max_prev = chain_image_down[k + nk * new_index];
						}
					}
					new_index = index + nx;
					if (bin_input_image[new_index] == 1 && chain_image_down[k + nk * new_index] > max_prev) {
						max_prev = chain_image_down[k + nk * new_index];
					}
					if (x < nx - 1) {
						new_index = index + nx + 1;
						if (bin_input_image[new_index] == 1 && chain_image_down[k + nk * new_index] > max_prev) {
							max_prev = chain_image_down[k + nk * new_index];
						}
					}

					/* Update chain length? */
					if (max_prev + 1 < chain_image_down[k + nk * index]) {
						// Update chain length
						chain_image_down[k + nk * index] = max_prev + 1;

						// Propagate changes to output
						if (bin_input_image[index]) {
							char new_bin_output_flag = 
								(chain_image_up[(K - k) + nk * index] + chain_image_down[k + nk * index] + 1 >= L);
							// Did we cross the threshold?
							if (bin_output_image_array[(K - k) + nk * index] && !new_bin_output_flag) {
								// Update flags
								bin_output_image_array[(K - k) + nk * index] = 0;
								--bin_output_image_count[index];
								// Did this extinguish the last path?
								if (bin_output_image_count[index] == 0) {
									// Write to output
									output_image[index] = threshold;
								}
							}
						} else {
							if (K - 1 - k >= 0) {
								char new_bin_output_flag = 
									(chain_image_up[(K - 1 - k) + nk * index] + chain_image_down[k + nk * index] + 1 >= L);
								// Did we cross the threshold?
								if (bin_output_image_array[(K - 1 - k) + nk * index] && !new_bin_output_flag) {
									// Update flags
									bin_output_image_array[(K - 1 - k) + nk * index] = 0;
									--bin_output_image_count[index];
									// Did this extinguish the last path?
									if (bin_output_image_count[index] == 0) {
										// Write to output
										output_image[index] = threshold;
									}
								}
							}
							// Else, this pixel has already been removed
						}

						/* Propagate changes by enqueueing upward neighbours */
						if (y > 0) {
							// Same layer
							if (x > 0) {
								new_index = index - nx - 1;
								if (!in_queue_up[k + nk * new_index]) {
									new_row_queue_cur_k.push_back(x - 1);
									in_queue_up[k + nk * new_index] = 1;
								}
							}
							new_index = index - nx;
							if (!in_queue_up[k + nk * new_index]) {
								new_row_queue_cur_k.push_back(x);
								in_queue_up[k + nk * new_index] = 1;
							}
							if (x < nx - 1) {
								new_index = index - nx + 1;
								if (!in_queue_up[k + nk * new_index]) {
									new_row_queue_cur_k.push_back(x + 1);
									in_queue_up[k + nk * new_index] = 1;
								}
							}
							// Down one layer
							if (k < K) {
								if (x > 0) {
									new_index = index - nx - 1;
									if (!in_queue_up[k + 1 + nk * new_index]) {
										new_row_queue_next_k.push_back(x - 1);
										in_queue_up[k + 1 + nk * new_index] = 1;
									}
								}
								new_index = index - nx;
								if (!in_queue_up[k + 1 + nk * new_index]) {
									new_row_queue_next_k.push_back(x);
									in_queue_up[k + 1 + nk * new_index] = 1;
								}
								if (x < nx - 1) {
									new_index = index - nx + 1;
									if (!in_queue_up[k + 1 + nk * new_index]) {
										new_row_queue_next_k.push_back(x + 1);
										in_queue_up[k + 1 + nk * new_index] = 1;
									}
								}
							}
						}
					}
				}
				// Wipe old queue
				row_queue.resize(0);

				/* Merge new row queues into existing queues */
				if (y - 1 >= 0) {
					if (new_row_queue_cur_k.size() > 0) {
						path_queue_up.merge_row(new_row_queue_cur_k, k, y - 1);
					}
					if (new_row_queue_next_k.size() > 0) {
						path_queue_up.merge_row(new_row_queue_next_k, k + 1, y - 1);
					}
				}
			}
		}
#ifdef DEBUGGING
		cout << "UPWARD SWEEP - after" << endl;
		dump_state(
			nx,									// Image dimensions etc.
			ny,
			nk,
			path_queue_up,						// Queueing structures
			path_queue_down,
			in_queue_up,
			in_queue_down,
			bin_input_image,					// Input/output binary images
			bin_output_image_array,
			bin_output_image_count,
			chain_image_up,						// Up/down chain lengths
			chain_image_down
		);
#endif
	}

	/* Free allocated memory */
	free((void *)in_queue_up);
	free((void *)in_queue_down);

	free((void *)bin_input_image);
	
	free((void *)chain_image_up);
	free((void *)chain_image_down);

	free((void *)bin_output_image_array);
	free((void *)bin_output_image_count);

	return 0;
}




/* A path opening in the ++ diagonal direction.
	Conjugate with flip to perform +- diagonal path openings
*/
static int diag_pathopen(
	PATHOPEN_PIX_TYPE * input_image,					/* The input image */
	int * sorted_indices,								/* Monotonic transform to [0, 1, ...] of input image */
	int nx, int ny,										/* Image dimensions */
	int L,												/* The threshold line length */
	int K,												/* The maximum gap number */
	PATHOPEN_PIX_TYPE * output_image					/* Output image */
)
{
	int i, k, x, y, index, new_index, sort_index, num_pixels;

	/************************************** Allocation **********************************************/
	num_pixels = nx * ny;
	int nk = K + 1;

	/* Construct queueing system */
	Path_Queue path_queue_up(nk, ny, nx);
	Path_Queue path_queue_down(nk, ny, nx);

	/* Dynamic binary input image */
	char * bin_input_image = (char *)malloc(num_pixels * sizeof(char));

	/* in_queue flags */
	char * in_queue_up = (char *)malloc(num_pixels * nk * sizeof(char));
	char * in_queue_down = (char *)malloc(num_pixels * nk * sizeof(char));

	/* Chain length images [k + nk * pixel_index].  These don't include the current pixel. */
	int * chain_image_up = (int *)malloc(num_pixels * nk * sizeof(int));
	int * chain_image_down = (int *)malloc(num_pixels * nk * sizeof(int));

	// At each pixel, we store the vector of binary outputs indexed by gap number of upward chain
	char * bin_output_image_array = (char *)malloc(num_pixels * nk * sizeof(char));
	// Also count the vector of binary outputs, to note when they are all extinguished (boolean PQ!)
	char * bin_output_image_count = (char *)malloc(num_pixels * sizeof(char));

	/************************************** Initialisation **********************************************/
	/* Dynamic binary threshold image is initially all 1's */
	memset(bin_input_image, 1, num_pixels * sizeof(char));

	/* Queue initially empty */
	memset(in_queue_up, 0, num_pixels * nk * sizeof(char));
	memset(in_queue_down, 0, num_pixels * nk * sizeof(char));

	/* Initialise the chain lengths */
	for (y = 0; y < ny; ++y) {
	for (x = 0; x < nx; ++x) {
		int index = x + nx * y;
		int up_length = x + y;
		int down_length = ((nx - 1) - x) + ((ny - 1) - y);

#ifdef PATHOPEN_LENGTH_HEURISTIC
		/* Threshold at L - 1 */
		if (up_length > L - 1) up_length = L - 1;
		if (down_length > L - 1) down_length = L - 1;
#endif

		for (k = 0; k < nk; ++k)
			chain_image_up[k + nk * index] = up_length;
		for (k = 0; k < nk; ++k)
			chain_image_down[k + nk * index] = down_length;
	}
	}

	/* Binary output vector at each pixel */
	memset(bin_output_image_array, 1, num_pixels * nk * sizeof(char));
	memset(bin_output_image_count, nk, num_pixels * sizeof(char));

	/* Set output image to default value (0 here) */
	memset(output_image, 0, num_pixels * sizeof(PATHOPEN_PIX_TYPE));
	/****************************************************************************************************/

	/* For each threshold from smallest to largest */
	sort_index = 0;
	while(sort_index < num_pixels) {
		PATHOPEN_PIX_TYPE threshold;

		/*********************************** Process threshold pixels *****************************************/
		threshold = input_image[sorted_indices[sort_index]];
#ifdef DEBUGGING
		cout << "Threshold = " << (int)threshold << endl;
#endif
		while(input_image[sorted_indices[sort_index]] == threshold) {
			/* Collect into rows for enqueueing */
			vector< vector<PIXEL_INDEX_TYPE> > new_row_queue_down(nk);
			vector< vector<PIXEL_INDEX_TYPE> > new_row_queue_right(nk);
			vector< vector<PIXEL_INDEX_TYPE> > new_row_queue_up(nk);
			vector< vector<PIXEL_INDEX_TYPE> > new_row_queue_left(nk);
			int row_y = sorted_indices[sort_index] / nx;
#ifdef DEBUGGING
			cout << "y = " << row_y << endl;
#endif

			while(
				input_image[sorted_indices[sort_index]] == threshold 
				&& 
				sorted_indices[sort_index] / nx == row_y
			) {
				/* Extract index and coordinates */
				index = sorted_indices[sort_index];
				y = index / nx;
				x = index % nx;

#ifdef DEBUGGING
				cout << "\tx = " << x << endl;
#endif

				/* Directly perform the changes to this threshold pixel */
				if (bin_input_image[index]) {
#ifdef DEBUGGING
					cout << "\tRemoving..." << endl;
					cout << "\tcount = " << (int)bin_output_image_count[index] << endl;
#endif
					// Remove this pixel from the binary input image
					bin_input_image[index] = 0;

					// Update the outputs
					if (bin_output_image_count[index] > 0) {
						// Update the output flags for each gap index (and count)
						bin_output_image_count[index] = 0;
						for (k = 0; k < K; ++k) {
							bin_output_image_array[k + nk * index] = 
								(chain_image_up[k + nk * index] + chain_image_down[(K - 1 - k) + nk * index] + 1) >= L;
							bin_output_image_count[index] += bin_output_image_array[k + nk * index];
						}

						// If all paths have been extinguished, update output
						if (bin_output_image_count[index] == 0) {
							output_image[index] = threshold;
						}
					}

					/* Enqueue downward pixels for update */
					if (y < ny - 1) {
						new_index = index + nx;
						// Enqueue this pixel for all k
						for (k = 0; k < nk; ++k) {
							if (!in_queue_down[k + nk * new_index]) {
								in_queue_down[k + nk * new_index] = 1;
								new_row_queue_down[k].push_back(x);
							}
						}
						if (x < nx - 1) {
							new_index = index + nx + 1;
							
							// Enqueue this pixel for all k
							for (k = 0; k < nk; ++k) {
								if (!in_queue_down[k + nk * new_index]) {
									in_queue_down[k + nk * new_index] = 1;
									new_row_queue_down[k].push_back(x + 1);
								}
							}
						}
					}
					if (x < nx - 1) {
						new_index = index + 1;
						// Enqueue this pixel for all k
						for (k = 0; k < nk; ++k) {
							if (!in_queue_down[k + nk * new_index]) {
								in_queue_down[k + nk * new_index] = 1;
								new_row_queue_right[k].push_back(x + 1);
							}
						}
					}

					/* Enqueue upward pixels for update */
					if (y > 0) {
						if (x > 0) {
							new_index = index - nx - 1;
							
							// Enqueue this pixel for all k
							for (k = 0; k < nk; ++k) {
								if (!in_queue_up[k + nk * new_index]) {
									in_queue_up[k + nk * new_index] = 1;
									new_row_queue_up[k].push_back(x - 1);
								}
							}
						}
						new_index = index - nx;
						// Enqueue this pixel for all k
						for (k = 0; k < nk; ++k) {
							if (!in_queue_up[k + nk * new_index]) {
								in_queue_up[k + nk * new_index] = 1;
								new_row_queue_up[k].push_back(x);
							}
						}
					}
					if (x > 0) {
						new_index = index - 1;
						
						// Enqueue this pixel for all k
						for (k = 0; k < nk; ++k) {
							if (!in_queue_up[k + nk * new_index]) {
								in_queue_up[k + nk * new_index] = 1;
								new_row_queue_left[k].push_back(x - 1);
							}
						}
					}
				}
				
				// Next pixel in row			
				++sort_index;
				if (sort_index >= num_pixels) break;
			}
#ifdef DEBUGGING
			cout << "Merging row queues" << endl;
#endif

			// Merge row queues into path queues
#ifdef DEBUGGING
			cout << "Down queues" << endl;
#endif
			// Down
			if (row_y + 1 < ny) {
				for (k = 0; k < nk; ++k) {
					if (new_row_queue_down[k].size() > 0) {
#ifdef DEBUGGING
						cout << "k = " << k << ":";
						for (i = 0; i < new_row_queue_down[k].size(); ++i) {
							cout << new_row_queue_down[k][i] << " ";
						}
						cout << endl;
#endif
						path_queue_down.merge_row(new_row_queue_down[k], k, row_y + 1);
					}
				}
			}
			// Right
			for (k = 0; k < nk; ++k) {
				if (new_row_queue_right[k].size() > 0) {
#ifdef DEBUGGING
				cout << "k = " << k << ":";
				for (i = 0; i < new_row_queue_right[k].size(); ++i) {
					cout << new_row_queue_right[k][i] << " ";
				}
				cout << endl;
#endif
					path_queue_down.merge_row(new_row_queue_right[k], k, row_y);
				}
			}


#ifdef DEBUGGING
			cout << "Up queues" << endl;
#endif
			// Up
			if (row_y - 1 >= 0) {
				for (k = 0; k < nk; ++k) {
					if (new_row_queue_up[k].size() > 0) {
#ifdef DEBUGGING
						cout << "k = " << k << ":";
						for (i = 0; i < new_row_queue_up[k].size(); ++i) {
							cout << new_row_queue_up[k][i] << " ";
						}
						cout << endl;
#endif
						path_queue_up.merge_row(new_row_queue_up[k], k, row_y - 1);
					}
				}
			}
			// Left
			for (k = 0; k < nk; ++k) {
				if (new_row_queue_left[k].size() > 0) {
#ifdef DEBUGGING
					cout << "k = " << k << ":";
					for (i = 0; i < new_row_queue_left[k].size(); ++i) {
						cout << new_row_queue_left[k][i] << " ";
					}
					cout << endl;
#endif
					path_queue_up.merge_row(new_row_queue_left[k], k, row_y);
				}
			}

#ifdef DEBUGGING
			cout << endl;
#endif
			if (sort_index >= num_pixels) break;
		}

		/*************************************** Downward sweep *********************************************/
		/* Propagate changes at current threshold down the image */
#ifdef DEBUGGING
		cout << "DOWNWARD SWEEP - before" << endl;
		dump_state(
			nx,									// Image dimensions etc.
			ny,
			nk,
			path_queue_up,						// Queueing structures
			path_queue_down,
			in_queue_up,
			in_queue_down,
			bin_input_image,					// Input/output binary images
			bin_output_image_array,
			bin_output_image_count,
			chain_image_up,						// Up/down chain lengths
			chain_image_down
		);
#endif

		for (k = 0; k < nk; ++k) {
#ifdef DEBUGGING
			cout << "k = " << k << endl;
#endif
			for (y = 0; y < ny; ++y) {
				vector<PIXEL_INDEX_TYPE> & row_queue = path_queue_down.q[k][y];
				if (row_queue.size() == 0) continue;


#ifdef DEBUGGING
				cout << "\ty = " << y << endl;
#endif

				bool right_queue = false;
				vector<PIXEL_INDEX_TYPE> cur_row_queue_next_k;
				vector<PIXEL_INDEX_TYPE> new_row_queue_cur_k;
				vector<PIXEL_INDEX_TYPE> new_row_queue_next_k;

				/* Perform updates on points in row_queue, propagating changes to the new row queues */
				i = 0;
				x = row_queue[i];
				while (true) {
					/* Extract x-coordinate and pixel index */
					index = x + nx * y;

#ifdef DEBUGGING
					cout << "\t\tx = " << x << endl;
#endif

					/* Unflag -> no longer in queue */
					in_queue_down[k + nk * index] = 0;

					/* Update chain length from upward neighbours */
					int max_prev = -1;
					// Previous level - accept a gap
					if (k > 0) {
						if (y > 0) {
							if (x > 0) {
								new_index = index - nx - 1;
								if (chain_image_up[k - 1 + nk * new_index] > max_prev) {
									max_prev = chain_image_up[k - 1 + nk * new_index];
								}
							}
							new_index = index - nx;
							if (chain_image_up[k - 1 + nk * new_index] > max_prev) {
								max_prev = chain_image_up[k - 1 + nk * new_index];
							}
						}
						if (x > 0) {
							new_index = index - 1;
							if (chain_image_up[k - 1 + nk * new_index] > max_prev) {
								max_prev = chain_image_up[k - 1 + nk * new_index];
							}
						}
					}
					// Current level - no gap allowed
					if (y > 0) {
						if (x > 0) {
							new_index = index - nx - 1;
							if (bin_input_image[new_index] == 1 && chain_image_up[k + nk * new_index] > max_prev) {
								max_prev = chain_image_up[k + nk * new_index];
							}
						}
						new_index = index - nx;
						if (bin_input_image[new_index] == 1 && chain_image_up[k + nk * new_index] > max_prev) {
							max_prev = chain_image_up[k + nk * new_index];
						}
					}
					if (x > 0) {
						new_index = index - 1;
						if (bin_input_image[new_index] == 1 && chain_image_up[k + nk * new_index] > max_prev) {
							max_prev = chain_image_up[k + nk * new_index];
						}
					}

					/* Update chain length? */
					if (max_prev + 1 < chain_image_up[k + nk * index]) {
#ifdef DEBUGGING
						cout << "New chain length is " << chain_image_up[k + nk * index] << endl;
#endif
						// Update chain length
						chain_image_up[k + nk * index] = max_prev + 1;

						// Propagate changes to output
						if (bin_input_image[index]) {
							char new_bin_output_flag = 
								(chain_image_up[k + nk * index] + chain_image_down[(K - k) + nk * index] + 1 >= L);
							// Did we cross the threshold?
							if (bin_output_image_array[k + nk * index] && !new_bin_output_flag) {
								// Clear the flag
								bin_output_image_array[k + nk * index] = 0;
								--bin_output_image_count[index];
								// Did this extinguish the last path?
								if (bin_output_image_count[index] == 0) {
									// Write to output
									output_image[index] = threshold;
								}
							}
						} else {
							if (K - 1 - k >= 0) {
								char new_bin_output_flag = 
									(chain_image_up[k + nk * index] + chain_image_down[(K - 1 - k) + nk * index] + 1 >= L);
								// Did we cross the threshold?
								if (bin_output_image_array[k + nk * index] && !new_bin_output_flag) {
									// Clear the flag
									bin_output_image_array[k + nk * index] = 0;
									--bin_output_image_count[index];
									// Did this extinguish the last path?
									if (bin_output_image_count[index] == 0) {
										// Write to output
										output_image[index] = threshold;
									}
								}
							}
							// Else, this pixel has already been removed
						}

						/* Propagate changes by enqueueing downward neighbours */
						if (y < ny - 1) {
							// Same layer
							new_index = index + nx;
							if (!in_queue_down[k + nk * new_index]) {
								new_row_queue_cur_k.push_back(x);
								in_queue_down[k + nk * new_index] = 1;
							}
							if (x < nx - 1) {
								new_index = index + nx + 1;
								if (!in_queue_down[k + nk * new_index]) {
									new_row_queue_cur_k.push_back(x + 1);
									in_queue_down[k + nk * new_index] = 1;
								}
							}
						}
						if (x < nx - 1) {
							new_index = index + 1;
							if (!in_queue_down[k + nk * new_index]) {
								right_queue = true;
								in_queue_down[k + nk * new_index] = 1;
							}
						}

						// Down one layer
						if (k < K) {
							if (y < ny - 1) {
								new_index = index + nx;
								if (!in_queue_down[k + 1 + nk * new_index]) {
									new_row_queue_next_k.push_back(x);
									in_queue_down[k + 1 + nk * new_index] = 1;
								}
								if (x < nx - 1) {
									new_index = index + nx + 1;
									if (!in_queue_down[k + 1 + nk * new_index]) {
										new_row_queue_next_k.push_back(x + 1);
										in_queue_down[k + 1 + nk * new_index] = 1;
									}
								}
							}
							if (x < nx - 1) {
								new_index = index + 1;
								if (!in_queue_down[k + 1 + nk * new_index]) {
									cur_row_queue_next_k.push_back(x + 1);
									in_queue_down[k + 1 + nk * new_index] = 1;
								}
							}
						}
					}

					/* Select next x */
					if (right_queue) {
						// Go across
						right_queue = false;
						++x;
						// Test halting condition
						if (x > nx - 1) break;

						// Consume the row queue?
						if (i + 1 < row_queue.size()) {
							if (row_queue[i + 1] == x)
								++i;
						}
					} else {
						// Halting condition
						if (i + 1 >= row_queue.size()) {
							break;
						} else {
							++i;
							x = row_queue[i];
						}
					}
				}
				// Wipe old queue
				row_queue.resize(0);

				/* Merge new row queues into existing queues */
				if (y + 1 < ny) {
					if (new_row_queue_cur_k.size() > 0) {
						path_queue_down.merge_row(new_row_queue_cur_k, k, y + 1);
					}
					if (new_row_queue_next_k.size() > 0) {
						path_queue_down.merge_row(new_row_queue_next_k, k + 1, y + 1);
					}
				}
				if (cur_row_queue_next_k.size() > 0) {
					path_queue_down.merge_row(cur_row_queue_next_k, k + 1, y);
				}
			}
		}
#ifdef DEBUGGING
		cout << "DOWNWARD SWEEP - after" << endl;
		dump_state(
			nx,									// Image dimensions etc.
			ny,
			nk,
			path_queue_up,						// Queueing structures
			path_queue_down,
			in_queue_up,
			in_queue_down,
			bin_input_image,					// Input/output binary images
			bin_output_image_array,
			bin_output_image_count,
			chain_image_up,						// Up/down chain lengths
			chain_image_down
		);
#endif

		/*************************************** Upward sweep *********************************************/
#ifdef DEBUGGING
		cout << "UPWARD SWEEP - before" << endl;
		dump_state(
			nx,									// Image dimensions etc.
			ny,
			nk,
			path_queue_up,						// Queueing structures
			path_queue_down,
			in_queue_up,
			in_queue_down,
			bin_input_image,					// Input/output binary images
			bin_output_image_array,
			bin_output_image_count,
			chain_image_up,						// Up/down chain lengths
			chain_image_down
		);
#endif

		/* Propagate changes at current threshold up the image */
		for (k = 0; k < nk; ++k) {
#ifdef DEBUGGING
			cout << "k = " << k << endl;
#endif
			for (y = ny - 1; y >= 0; --y) {
				vector<PIXEL_INDEX_TYPE> & row_queue = path_queue_up.q[k][y];
				if (row_queue.size() == 0) continue;

#ifdef DEBUGGING
				cout << "\ty = " << y << endl;
#endif

				bool left_queue = false;
				vector<PIXEL_INDEX_TYPE> cur_row_queue_next_k;
				vector<PIXEL_INDEX_TYPE> new_row_queue_cur_k;
				vector<PIXEL_INDEX_TYPE> new_row_queue_next_k;

				/* Perform updates on points in row_queue, propagating changes to the new row queues */
				i = row_queue.size() - 1;
				x = row_queue[i];
				while(true) {
					/* Extract x-coordinate and pixel index */
					index = x + nx * y;

#ifdef DEBUGGING
					cout << "\t\tx = " << x << endl;
#endif

					/* Unflag -> no longer in queue */
					in_queue_up[k + nk * index] = 0;

					/* Update chain length from downward neighbours */
					// Note: Only y < ny - 1 may be 'updated', so we are assured of the existence of previous neighbours!
					int max_prev = -1;
					// Previous level - accept a gap
					if (k > 0) {
						if (y < ny - 1) {
							new_index = index + nx;
							if (chain_image_down[k - 1 + nk * new_index] > max_prev) {
								max_prev = chain_image_down[k - 1 + nk * new_index];
							}
							if (x < nx - 1) {
								new_index = index + nx + 1;
								if (chain_image_down[k - 1 + nk * new_index] > max_prev) {
									max_prev = chain_image_down[k - 1 + nk * new_index];
								}
							}
						}
						if (x < nx - 1) {
							new_index = index + 1;
							if (chain_image_down[k - 1 + nk * new_index] > max_prev) {
								max_prev = chain_image_down[k - 1 + nk * new_index];
							}
						}
					}
					// Current level - no gap allowed
					if (y < ny - 1) {
						new_index = index + nx;
						if (bin_input_image[new_index] == 1 && chain_image_down[k + nk * new_index] > max_prev) {
							max_prev = chain_image_down[k + nk * new_index];
						}
						if (x < nx - 1) {
							new_index = index + nx + 1;
							if (bin_input_image[new_index] == 1 && chain_image_down[k + nk * new_index] > max_prev) {
								max_prev = chain_image_down[k + nk * new_index];
							}
						}
					}
					if (x < nx - 1) {
						new_index = index + 1;
						if (bin_input_image[new_index] == 1 && chain_image_down[k + nk * new_index] > max_prev) {
							max_prev = chain_image_down[k + nk * new_index];
						}
					}

					/* Update chain length? */
					if (max_prev + 1 < chain_image_down[k + nk * index]) {
						// Update chain length
						chain_image_down[k + nk * index] = max_prev + 1;

						// Propagate changes to output
						if (bin_input_image[index]) {
							char new_bin_output_flag = 
								(chain_image_up[(K - k) + nk * index] + chain_image_down[k + nk * index] + 1 >= L);
							// Did we cross the threshold?
							if (bin_output_image_array[(K - k) + nk * index] && !new_bin_output_flag) {
								// Update flags
								bin_output_image_array[(K - k) + nk * index] = 0;
								--bin_output_image_count[index];
								// Did this extinguish the last path?
								if (bin_output_image_count[index] == 0) {
									// Write to output
									output_image[index] = threshold;
								}
							}
						} else {
							if (K - 1 - k >= 0) {
								char new_bin_output_flag = 
									(chain_image_up[(K - 1 - k) + nk * index] + chain_image_down[k + nk * index] + 1 >= L);
								// Did we cross the threshold?
								if (bin_output_image_array[(K - 1 - k) + nk * index] && !new_bin_output_flag) {
									// Update flags
									bin_output_image_array[(K - 1 - k) + nk * index] = 0;
									--bin_output_image_count[index];
									// Did this extinguish the last path?
									if (bin_output_image_count[index] == 0) {
										// Write to output
										output_image[index] = threshold;
									}
								}
							}
							// Else, this pixel has already been removed
						}

						/* Propagate changes by enqueueing upward neighbours */
						if (y > 0) {
							// Same layer
							if (x > 0) {
								new_index = index - nx - 1;
								if (!in_queue_up[k + nk * new_index]) {
									new_row_queue_cur_k.push_back(x - 1);
									in_queue_up[k + nk * new_index] = 1;
								}
							}
							new_index = index - nx;
							if (!in_queue_up[k + nk * new_index]) {
								new_row_queue_cur_k.push_back(x);
								in_queue_up[k + nk * new_index] = 1;
							}
						}
						if (x > 0) {
							new_index = index - 1;
							if (!in_queue_up[k + nk * new_index]) {
								left_queue = true;
								in_queue_up[k + nk * new_index] = 1;
							}
						}

						// Down one layer
						if (k < K) {
							if (y > 0) {
								if (x > 0) {
									new_index = index - nx - 1;
									if (!in_queue_up[k + 1 + nk * new_index]) {
										new_row_queue_next_k.push_back(x - 1);
										in_queue_up[k + 1 + nk * new_index] = 1;
									}
								}
								new_index = index - nx;
								if (!in_queue_up[k + 1 + nk * new_index]) {
									new_row_queue_next_k.push_back(x);
									in_queue_up[k + 1 + nk * new_index] = 1;
								}
							}
							if (x > 0) {
								new_index = index - 1;
								if (!in_queue_up[k + 1 + nk * new_index]) {
									cur_row_queue_next_k.push_back(x - 1);
									in_queue_up[k + 1 + nk * new_index] = 1;
								}
							}
						}
					}

					/* Select next x */
					if (left_queue) {
						// Go across
						left_queue = false;
						--x;
						// Test halting condition
						if (x < 0) break;

						// Consume the row queue?
						if (i - 1 >= 0) {
							if (row_queue[i - 1] == x)
								--i;
						}
					} else {
						// Halting condition
						if (i - 1 < 0) {
							break;
						} else {
							--i;
							x = row_queue[i];
						}
					}
				}
				// Wipe old queue
				row_queue.resize(0);

				/* Merge new row queues into existing queues */
				if (y - 1 >= 0) {
					if (new_row_queue_cur_k.size() > 0) {
						reverse(new_row_queue_cur_k.begin(), new_row_queue_cur_k.end());
						path_queue_up.merge_row(new_row_queue_cur_k, k, y - 1);
					}
					if (new_row_queue_next_k.size() > 0) {
						reverse(new_row_queue_next_k.begin(), new_row_queue_next_k.end());
						path_queue_up.merge_row(new_row_queue_next_k, k + 1, y - 1);
					}
				}
				if (cur_row_queue_next_k.size() > 0) {
					reverse(cur_row_queue_next_k.begin(), cur_row_queue_next_k.end());
					path_queue_up.merge_row(cur_row_queue_next_k, k + 1, y);
				}
			}
		}
#ifdef DEBUGGING
		cout << "UPWARD SWEEP - after" << endl;
		dump_state(
			nx,									// Image dimensions etc.
			ny,
			nk,
			path_queue_up,						// Queueing structures
			path_queue_down,
			in_queue_up,
			in_queue_down,
			bin_input_image,					// Input/output binary images
			bin_output_image_array,
			bin_output_image_count,
			chain_image_up,						// Up/down chain lengths
			chain_image_down
		);
#endif
	}

	/* Free allocated memory */
	free((void *)in_queue_up);
	free((void *)in_queue_down);

	free((void *)bin_input_image);
	
	free((void *)chain_image_up);
	free((void *)chain_image_down);

	free((void *)bin_output_image_array);
	free((void *)bin_output_image_count);

	return 0;
}


/* dump_state:
	Nifty debugging function 
	- dump the state of the algorithm to the terminal for verification.
*/
void dump_state(
	int nx,									// Image dimensions etc.
	int ny,
	int nk,
	Path_Queue & path_queue_up,				// Queueing structures
	Path_Queue & path_queue_down,
	char * in_queue_up,
	char * in_queue_down,
	char * bin_input_image,					// Input/output binary images
	char * bin_output_image_array,
	char * bin_output_image_count,
	int * chain_image_up,					// Up/down chain lengths
	int * chain_image_down
)
{
	int i, x, y, k;
	int K = nk - 1;

	cout << endl << "STATE DUMP" << endl;

	// Observe the input binary image
	cout << "bin_input_image" << endl;
	for (y = 0; y < ny; ++y) {
		cout << "\t" << y << ":";
		for (x = 0; x < nx; ++x) {
			cout << (int)bin_input_image[x + nx * y] << " ";
		}
		cout << endl;
	}
	// Observe the bin_output_image_array
	cout << "bin_output_image_array" << endl;
	for (k = 0; k < nk; ++k) {
		cout << "k = " << k << endl;
		for (y = 0; y < ny; ++y) {
			cout << "\t" << y << ":";
			for (x = 0; x < nx; ++x) {
				cout << (int)bin_output_image_array[k + nk * (x + nx * y)] << " ";
			}
			cout << endl;
		}
	}
	// Observe the bin_output_image_count
	cout << "bin_output_image_count" << endl;
	for (y = 0; y < ny; ++y) {
		cout << "\t" << y << ":";
		for (x = 0; x < nx; ++x) {
			cout << (int)bin_output_image_count[x + nx * y] << " ";
		}
		cout << endl;
	}

	// Observe the chain images
	// Down
	cout << "chain_image_down:" << endl;
	for (k = 0; k < nk; ++k) {
		cout << "k = " << k << endl;
		for (y = 0; y < ny; ++y) {
			cout << "\t" << y << ":";
			for (x = 0; x < nx; ++x) {
				cout << chain_image_down[k + nk * (x + nx * y)] << " ";
			}
			cout << endl;
		}
	}
	// Up
	cout << "chain_image_up:" << endl;
	for (k = 0; k < nk; ++k) {
		cout << "k = " << k << endl;
		for (y = 0; y < ny; ++y) {
			cout << "\t" << y << ":";
			for (x = 0; x < nx; ++x) {
				cout << chain_image_up[k + nk * (x + nx * y)] << " ";
			}
			cout << endl;
		}
	}
		
	// Observe the path queues
	/*
	cout << "path_queue_down:" << endl;
	for (k = 0; k < nk; ++k) {
		cout << "k = " << k << endl;
		for (y = 0; y < ny; ++y) {
			vector<PIXEL_INDEX_TYPE> & q = path_queue_down.q[k][y];
			if (q.size() == 0) continue;

			cout << "\t" << y << ":";
			for (i = 0; i < q.size(); ++i) {
				cout <<	q[i] << " ";
			}
			cout << endl;
		}
	}
	*/
	cout << "in_queue_down" << endl;
	for (k = 0; k < nk; ++k) {
		cout << "k = " << k << endl;
		for (y = 0; y < ny; ++y) {
			cout << "\t" << y << ":";
			for (x = 0; x < nx; ++x) {
				cout << (int)in_queue_down[k + nk * (x + nx * y)] << " ";
			}
			cout << endl;
		}
	}
	/*
	cout << "path_queue_up:" << endl;
	for (k = 0; k < nk; ++k) {
		cout << "k = " << k << endl;
		for (y = 0; y < ny; ++y) {
			vector<PIXEL_INDEX_TYPE> & q = path_queue_up.q[k][y];
			if (q.size() == 0) continue;

			cout << "\t" << y << ":";
			for (i = 0; i < q.size(); ++i) {
				cout <<	q[i] << " ";
			}
			cout << endl;
		}
	}
	*/
	cout << "in_queue_up" << endl;
	for (k = 0; k < nk; ++k) {
		cout << "k = " << k << endl;
		for (y = 0; y < ny; ++y) {
			cout << "\t" << y << ":";
			for (x = 0; x < nx; ++x) {
				cout << (int)in_queue_up[k + nk * (x + nx * y)] << " ";
			}
			cout << endl;
		}
	}

	/*
	cout << "getch()" << endl;
	getch();
	*/
}
