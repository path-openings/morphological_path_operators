/*
 *		File:		path_queue.h
 *
 *		Purpose:	Specialised queues for path openings
 *
 *		Author:		Ben Appleton
 *		Date:		30/06/2005
 *
 *              Latest modifications : Hugues Talbot	30 Nov 2009

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

 

/*		Notes:
 *			- Does not account for duplicates (handle outside)
 *			- Expects and maintains ascending order
 */

#ifndef PATH_QUEUE_H
#define PATH_QUEUE_H

#include <vector>
using namespace std;

// The type used to index pixels in general
typedef int PIXEL_INDEX_TYPE;
static const int PIXEL_INDEX_TYPE_MIN = 0;
static const int PIXEL_INDEX_TYPE_MAX = 0x7fffffff;

// The type used to store pixel data
typedef unsigned char PIXEL_TYPE;
static const int PIXEL_TYPE_MIN = 0;
static const int PIXEL_TYPE_MAX = 255;

class Path_Queue {
public:
	/* Data */
	// Queue data
	int num_gaps;
	int num_rows;
	int row_max_length;
	vector< vector< vector< PIXEL_INDEX_TYPE > > > q;

	/* Methods */
	Path_Queue(
		int num_gaps,
		int num_rows,
		int row_max_length
	);

	~Path_Queue();

	/* merge_row:
		Merge the given list of indices with the queue for the specified row.
		NOTE: This is the only means of insertion; 
			this is to emphasise that insertion of single elements is slow!

	*/
	void merge_row(
		vector<PIXEL_INDEX_TYPE> row,
		int k,
		int r
	);

	/* print_state:
		Debugging function - dumps the internal state
	*/
	void print_state();
};

#endif
