/*+++++++++++++++++ class VertexLists +++++++++++++++++++ */
#ifndef __vertexlists_h__
#define __vertexlists_h__

#include <iostream>
#include <vector>

#include <Exceptions/vertexlistsexception.h>

/*
 Implements a set of lists where each list contains
 vertices. It is implemented as an array. Each element
 contains the index of its precursor and successor in
 the array. Thus, a list is a double-linked way
 through the array. The heads of the list are saved
 in a head-list. Elements can be accessed directly
 via their index or by iterating through a list.
 Each element contains its acutual list where it belongs
 to.

 Elements can be deleted from lists, copied to other lists
 and lists can be merged as well.
*/

class VertexLists {
public:

	/* #### Type definitions #### */

	// type to save a index, we use short int size vertex indices are smaller than 65535
	typedef unsigned short pos_int_type;

	// each element in the matrix contains its precursor, successor and the lists it is in
	struct vertex_type
	{
		pos_int_type next;
		pos_int_type prev;
		pos_int_type list;
	};

	// a head element saves the number of elements in the list and the first element of it
	struct head_type
	{
		pos_int_type head;
		int size;
	};
	
	// vertex set type, to return a list of vertices
	typedef std::vector< pos_int_type > vertex_set_type;

	// list of elements in the array
	typedef std::vector< vertex_type > vertex_list_type;

	// head-list
	typedef std::vector< head_type > head_list_type;

	// function pointer to short IO function
	
	// reference to nothing, if element has e.g. no successor
	const pos_int_type nil;



	/* #### Constructor #### */

	// empty constructor
	VertexLists();

	// constructor which saves number of possible elements and number of possible lists
	VertexLists(int size, pos_int_type number_of_lists);

	/* ####  Destructor  #### */

	~VertexLists() = default;



	/* ####  public methods  #### */

	// returns number of possible elements
	int getSize();

	// returns number of lists	
	int getNumberOfLists();

	// return list number of element i
	pos_int_type getListNr(pos_int_type i);

	// return size of list "list"
	int getListSize(pos_int_type list);

	// return top element of list "list"
	pos_int_type getTop(pos_int_type list);

	// return all elements in list "list" in a vertex_set_type list
	vertex_set_type getList(pos_int_type list);

	// add element i to list "list", just possible if i is not member of any list
	void addToList(pos_int_type, pos_int_type list);

	// deletes element i from its list, just possible if i is member of any list
	void deleteFromList(pos_int_type i);

	// copy vertex i to list "list", just possible i is member of a list
	void copyToList(pos_int_type i, pos_int_type list);

	// merges two lists i,j, such that all elements from j are appended
	// to list to i and list j is empty at its end
	void mergeLists(pos_int_type i, pos_int_type j);


	/* ####  programmer's help functions #### */

		
	// print all lists on std output
	void printLists();

	// print list "list" on std output
	void printList(pos_int_type list);

private:
	/* ####  member variables  #### */
	// list of elements
	vertex_list_type _vertex_list;

	// list of heads (head pointers)
	head_list_type _heads;
	
	// number of lists
	pos_int_type _number_of_lists;

	// number of possible elements
	int _size;

};

#endif
