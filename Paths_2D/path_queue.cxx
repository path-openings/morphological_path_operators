/*
 *		File:		path_queue.cxx
 *
 *		Purpose:	Specialised queues for path openings
 *
 *		Author:		Ben Appleton
 *

  Copyright Benjamin Appleton and Hugues Talbot, July 2005

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

#pragma warning(disable : 4786) // Disable the silly 'warning' about templated variable names becoming too long

#include "path_queue.h"

#include <iostream>
using namespace std;

Path_Queue::Path_Queue(
	int num_gaps,
	int num_rows,
	int row_max_length
)
{
	int i;

	this->num_gaps = num_gaps;
	this->num_rows = num_rows;
	this->row_max_length = row_max_length;

	/* Set up the size of the path queue */
	q.resize(num_gaps);
	for (i = 0; i < num_gaps; ++i) {
		q[i].resize(num_rows);
	}
	// Individual rows are a vector list implementation; nothing to do
}

Path_Queue::~Path_Queue()
{
	// Nothing to do
}

/* merge_row:
	Merge the given list of indices with the queue for the specified row.
	Maintains ascending order.
*/
void Path_Queue::merge_row(
	vector<PIXEL_INDEX_TYPE> row,		// Assumed non-empty for efficiency
	int k,								// Gap number
	int r								// Row number
)
{
	vector<PIXEL_INDEX_TYPE> & old_row = this->q[k][r];

	/*  Shortcut */
	if (old_row.size() == 0) {
		old_row = row;
		return;
	}

	vector<PIXEL_INDEX_TYPE> new_row;		// Initially empty

	/* DEBUGGING */
	/*
	this->print_state();
	int j;
	for (j = 0; j < row.size(); ++j) {
		cout << row[j] << ", ";
	}
	cout << endl;
	*/
	/*************/

	/* Merge row and old_row into new_row
		- consumes them in competition
		- uses sentinels to detect end of row
	*/
	int old_row_index;
	int row_index;
	int old_row_size = old_row.size();
	int row_size = row.size();
	// Append sentinels
	row.push_back(PIXEL_INDEX_TYPE_MAX);
	old_row.push_back(PIXEL_INDEX_TYPE_MAX);
	for (
		old_row_index = 0, row_index = 0; 
		old_row_index < old_row_size || row_index < row_size;		// Halt when both rows are consumed

	) {
		if (old_row[old_row_index] < row[row_index]) {
			new_row.push_back(old_row[old_row_index]);
			++old_row_index;
		} else {
			new_row.push_back(row[row_index]);
			++row_index;
		}
	}
	// Remove the sentinel from 'row'
	row.pop_back();

	// Replace the old row by the new row
	old_row = new_row;
}

/* print_state:
	Debugging function - dumps the internal state
*/
void Path_Queue::print_state() {
	int k, r, i;
	
	cout << "print_state:" << endl;
	for (k = 0; k < num_gaps; ++k) {
		for (r = 0; r < num_rows; ++r) {
			vector<PIXEL_INDEX_TYPE> & row = q[k][r];
			int length = row.size();
			if (length > 0) {
				cout << "(" << k << ", " << r << "): ";
				for (i = 0; i < length; ++i) {
					cout << row[i];
					if (i < length - 1) cout << ", ";
				}
				cout << endl;
			}
		}
	}
}
