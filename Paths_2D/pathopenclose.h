/*
 * File:		pathopenclose.h
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
 pathopenclose.h
 ------

  DESCRIPTION:
  Path openings and closings for grayscale images.

  REFERENCE:
  Previous work: Hugues Talbot, Michael Buckley, Henk Heijmans (in no particular order)

  HISTORY:
  Created by Ben Appleton (July 2004)
  Contact: appleton@itee.uq.edu.au
**********************************************************************************************/

#ifndef PATHOPENCLOSE_H
#define PATHOPENCLOSE_H

extern "C" {
	#include "path_support.h"
}

/* Define DEBUG to enable printfs for debugging */
/* #define PATHOPEN_DEBUG */
/* #define PATHOPEN_DIAG_DEBUG */

int pathopen(
	PATHOPEN_PIX_TYPE * input_image,				/* The input image */
	int nx, int ny,									/* Image dimensions */
	int L,											/* The threshold line length */
	int K,											/* The maximum number of gaps in the path */
	PATHOPEN_PIX_TYPE * output_image				/* Output image */
);

#endif // PATHOPENCLOSE_H
