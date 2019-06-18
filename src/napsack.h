/*+++++++++++++++++ class Napsack +++++++++++++++++++ */
#ifndef __napsack_h__
#define __napsack_h__

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <array.h>

using std::vector;
using std::cout;
using std::endl;

/* This static class implements a Napsack like algorithm
   which works by dynamic programming and finds the optimal
   way to put tuples in the three buckets. For details see thesis.
   
   The calculation is aborted if the maximum already reaches a
   upper bound.
*/

class Napsack {

public:

	Napsack();	
	
	~Napsack();	
	
	/* ##### Type definitions ##### */
	
	// tuple (x,y) with up rounded and down rounded x and y
	struct Tupel {
		long long x_up;
		long long x_down;
		long long y_up; 
		long long y_down;
	};
	
	// field type in linear array, bool if it was used, max_y the maximal y entry in field x
	struct field_type {
		bool set;
		long long max_y;
	};
	
	
	// calculation of min-max napsack from given tuple list
	inline static long double getNapsackSolution(vector<Tupel> & Tupel_list, long long &matrix_x_size, long long &matrix_y_size, long double &delta_u, long double &delta_v, long double upper_bound, long long infinity)
	{	
		// create linear array
		Array<field_type> array = Array<field_type>(matrix_x_size);
		// create inf field
		field_type inf_field;

		// compute inital instance of array, just 0 is set to true
		compute_D_null(array, matrix_x_size, inf_field);

		// compute all following instance by adding the tuples to the buckets
		for (int i = 0; i < Tupel_list.size(); i++){
			// compute next instance by adding the next tuple
			long long max_value = compute_second_D(array, inf_field, Tupel_list[i].x_up, Tupel_list[i].x_down, Tupel_list[i].y_up, Tupel_list[i].y_down, matrix_x_size, infinity);
			
			// check if upper bound is already achieved, if so stop
			if (max_value+std::min(delta_u, delta_v) > upper_bound) {
				return max_value+std::min(delta_u,delta_v);
			}
		}
				
		// calculate min max of the last instance of the array
		long double max_value=-1;			
		double value;
		for (long long i = 0; i < array.size(); i++) {
			if (array[i].set) {
				value = std::min(compute_real_coordinate(i, matrix_x_size)+delta_u, array[i].max_y+delta_v);
				max_value = (value > max_value) ? value : max_value;
			}
		}
		max_value = ( inf_field.max_y + delta_v > max_value) ? (inf_field.max_y + delta_v)  : max_value;
		
		return max_value;
	}	


private:

	/* #### help functions #### */

	// calculate real coordinates
	inline static long long compute_real_coordinate (long long x, long long size)
	{
		return (x - (size/2));
	}


	// calculate matrix indices from real indices
	inline static long long compute_matrix_index (long long y, long long size)
	{
		return (y + (size/2));
	}


	// fill inital instance of array..by setting just 0 to true
	inline static void compute_D_null(Array<field_type> &array, long long x_size, field_type &inf_field)
	{
		// fill all field with 0, false
		field_type zero = {false, 0};
		for (long long i = 0; i < array.size(); i++) {
			array[i] = zero;
		}
		inf_field = zero;
		
		// finally set middle field...field 0 to true
		field_type one = {true, 0};
		array[compute_matrix_index(0,x_size)] = one;
	}
	

	// calculate following instance of array by adding the next tuple
	// also checking the maximal entry on the fly
	inline static long long compute_second_D(Array<field_type> &array, field_type &inf_field, long long x_up, long long x_down, long long y_up, long long y_down, long long x_size, long long infinity)
	{
		// intiate maximal entry with zero
		long long max_entry = 0;
		
		// create new array instance from old one...old positions which have been true are still true
		Array<field_type> new_array = array;
		long long min_v;

		// if x is -infinity then update just the inf field
		if (x_down < -infinity) {
			//std::cout << "new tuple is -inf with y=" << y_up << std::endl;
			// iterate over old array instance
			for (long long  i = 0; i < array.size(); i++) {
				// if old instance was true set the new instance at two positions to true
				if (array[i].set) {
					if (inf_field.set == false || array[i].max_y + y_up > inf_field.max_y) {
						field_type new_value = {true, array[i].max_y+ y_up};
						//std::cout << "inf field is " << array[i].max_y+ y_up  << std::endl;
						inf_field = new_value;
						// calculate min value of both, and check it this one is greater as maximum
						min_v =  array[i].max_y + y_up;
						max_entry = ( min_v > max_entry) ? min_v : max_entry;
					}
				}
			}
		} else {
			//std::cout << "new tuple x=" << x_down << " y=" << y_down << std::endl;
			// update field dependend on x-values

			// its always true that inf_field.max_y + y > inf_field.max_y
			//std::cout << "at first field inf is " << inf_field.max_y + y_up << std::endl;
			field_type new_value = {true, inf_field.max_y + y_up};
			inf_field = new_value;
			min_v = inf_field.max_y;
			max_entry = ( min_v > max_entry) ? min_v : max_entry;

			// iterate over old array instance
			for (long long  i = 0; i < array.size(); i++) {
				// if old instance was true set the new instance at two positions to true
				if (array[i].set) {
					// standard case
					if (array[i-x_down].set == false || array[i].max_y + y_up > array[i-x_down].max_y) {
						//std::cout << "field " << i-x_down << " is " << array[i].max_y+ y_up << std::endl;
						field_type new_value = {true, array[i].max_y+ y_up};
						new_array[(i-x_down)] = new_value;
						// calculate min value of both, and check it this one is greater as maximum
						min_v = std::min( compute_real_coordinate(i-x_down, x_size), array[i].max_y+ y_up);
						max_entry = ( min_v > max_entry) ? min_v : max_entry;
					}
					if (array[i+x_up].set == false || array[i].max_y - y_down > array[i+x_up].max_y) {
						//std::cout << "field " << i+x_up << " is " << array[i].max_y - y_down << std::endl;
						field_type new_value = {true, array[i].max_y - y_down};
						new_array[i+x_up] = new_value;
						// calculate min value of both, and check it this one is greater as maximum
						min_v = std::min( compute_real_coordinate(i+x_up, x_size), array[i].max_y - y_down);
						max_entry = ( min_v > max_entry) ? min_v : max_entry;
					}
					if (inf_field.set == false  || array[i].max_y + y_up > inf_field.max_y) {
						//std::cout << "field inf  is " << array[i].max_y + y_up << std::endl;
						field_type new_value = {true, array[i].max_y + y_up};
						inf_field = new_value;
						min_v =  array[i].max_y + y_up;
						max_entry = ( min_v > max_entry) ? min_v : max_entry;
					}

				}
			}
			
		}
				
		// save new array instance
		array = new_array;
		
		// and return maximum so far
		return max_entry;
	}

};
#endif
