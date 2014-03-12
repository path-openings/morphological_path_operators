/*
 * File:		pathopen.h
 *
 * Written by:		Ben Appleton
 *
 * Date:		February 2004
 
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
 pathopen.h
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "pathopenclose.h"
#include "path_queue.h"

extern "C" {
	#include "path_support.h"
}

#define PATHOPEN_LENGTH_HEURISTIC

/************************************* FUNCTION PROTOTYPES **************************************/
/* A path opening along the vertical direction.
Conjugate with transpose to perform horizontal path openings */
static int vert_pathopen(
	PATHOPEN_PIX_TYPE * input_image,					/* The input image */
	int * sorted_indices,								/* Monotonic transform to [0, 1, ...] of input image */
	int nx, int ny,										/* Image dimensions */
	int L,												/* The threshold line length */
	int K,												/* The maximum gap number */
	PATHOPEN_PIX_TYPE * output_image					/* Output image */
);

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
);

/* Debugging */
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
);
