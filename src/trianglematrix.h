/*+++++++++++++++++ class TriangleMatrix +++++++++++++++++++ */
#ifndef __trianglematrix_h__
#define __trianglematrix_h__

#include <iostream>

#include <array.h>

/*
Implements a generic triangle matrix as it is often used for
an undirected graph adjacency matrix.
It is possible to delete certain indices as if you would
use std::vector. However, a delted vertex is not deleted
for real. Instead a translation table for vertex indices is
used and just here an index is delted.

It is implemented using an Array object. To redirect in the 
right rows and columns it uses a offset matrix, which contains
in field i the number of elements before row i start.

Access is done via pos function and dirpos function.
Dirpos function is a direct access to an element, where the
first index always has to be greater than the second one.
Pos function gives access for two indices.

Note that both access function won't check whether you stay
in the triangle matrix or not. It also won't check for
self indication (meaning index1 = index2)!!
*/
template<class Type>
class TriangleMatrix {
public:

    /* ####  Constructors  #### */

    // empty constructor
    inline TriangleMatrix() : _size(0), _matrix(0) {
    }

    // initates matrix of given size with given inital values
    inline TriangleMatrix(int size, Type init_value) : _size(size) {
        // create matrix
        _matrix = Array<Type>(static_cast<int>((static_cast<double>(size) / 2.0) * (size - 1)), init_value);

        // create offset vector and index translation table
        _offset_array = Array<int>(size);
        _index_list = Array<int>(size);

        // fill offset vector and index translation table
        for (int i = 0; i < size; i++) {
            _index_list[i] = i;
            _offset_array[i] = static_cast<int>((static_cast<double>(i) / 2.0) * (i - 1.0));
        }
    }

    // copy constructor
    inline TriangleMatrix<Type>(TriangleMatrix<Type> const &new_matrix) : _size(new_matrix._size) {
        // create matrix
        _matrix = new_matrix._matrix;

        // create offset vector
        _offset_array = new_matrix._offset_array;
        _index_list = new_matrix._index_list;
    }


    // destructor
    inline ~TriangleMatrix() {
        _index_list.~Array<int>();
        _matrix.~Array<Type>();
        _offset_array.~Array<int>();
    }


    // return size of matrix
    inline int size() const {
        return _size;
    }


    // assignment function
    inline TriangleMatrix<Type> &operator=(const TriangleMatrix<Type> &right) {
        // assign all members
        _matrix = right._matrix;
        _offset_array = right._offset_array;
        _size = right._size;
        _index_list = right._index_list;

        return *this;
    }


    // return interal index for outerworld index i
    inline int getInternalIndex(int i) const {
        return _index_list[i];
    }

    // access function to access position i,k (i,k can be swapped)
    inline Type &pos(int i, int k) const {
        // return interal true index of i,k
        i = getInternalIndex(i);
        k = getInternalIndex(k);

        // swap i,k if necessary
        if (i < k) {
            int help = i;
            i = k;
            k = help;
        }

        // calculate index in internal array
        int index = _offset_array[i] + k;

        return _matrix[index];
    }


    // access function which gets position directly, condition i > k
    inline Type &dirpos(int i, int k) const {
        // return internal true index of i,k
        i = getInternalIndex(i);;
        k = getInternalIndex(k);;

        // calculate index in array and return it
        return _matrix[_offset_array[i] + k];
    }


    // deltes a vertex of the triangle matrix
    // note that vertex indices after index "index" are decreased by one
    // for the outside world
    inline void deleteIndex(int index) {
        // decrease virtually the index of all following indices
        for (int k = index; k < _index_list.size() - 1; k++) {
            _index_list[k] = _index_list[k + 1];
        }

        // decrease size variable of object
        _size--;
    }


private:
    /* ####  member variables  #### */

    // size of triangle matrix
    int _size;

    // index translation table
    Array<int> _index_list;

    // the actual triangle matrix uses an array
    Array<Type> _matrix;

    // offset list to saven row offset
    Array<int> _offset_array;


};

#endif
