/*+++++++++++++++++ class Array +++++++++++++++++++ */
#ifndef __array_h__
#define __array_h__

#include <iostream>

/*
Generic class for a array of static length.
Can be used instead of vector, since vector of std is slower
and larger in memory space.
It can be used equivalant to a normal array with the [] operator.

Note, as a normal array it won't check whether you stay
in bounds or not!!
*/

template <class Type> class Array {
public:

	/* ####  Constructors  #### */
	// empty constructor
	inline Array() : _size(0)
	{
		_array = NULL;
	}


	// initiate array with given size
	inline Array(int size) : _size(size)
	{
		// create matrix
		_array = new Type[size];
	}


	// initiate array with given size and given initial value
	inline Array(int size, Type init_value) : _size(size)
	{
		// create matrix
		_array = new Type[size];

		for (int i = 0; i < size; i++) {
		     _array[i] = init_value;
		}
	}


	// copy constructor
	inline Array<Type>(Array<Type> const &new_array) : _size(new_array._size)
	{
		_array = new Type[_size];

		for (int i = 0; i < _size; i++) {
		     _array[i] = new_array[i];
		}
	}


	// destructor
	inline ~Array()
	{
		delete[] _array;
		_array = NULL;
	}


	// return size of array
	inline int size() const
	{
		return _size;
	}


	// assignment of an array	
	inline Array<Type>& operator=(const Array<Type>& right)
	{
		// if size is not equal, a new array has to be created
		if (_size != right._size) {
		    // delete old array
		    delete[] _array;
		    _array = NULL;
		    // create a new array of necessary size
		    _size = right._size;
		    _array = new Type[_size];
		}

		// copy values
		for (int i = 0; i < _size; i++) {
		     _array[i] = right[i];
		}

		return *this;
	}


	// access operator
	inline Type& operator[](int index) const
	{
		return _array[index];
	}


private:
	/* ####  member variables  #### */
	int _size;
	Type *_array;
};

#endif
