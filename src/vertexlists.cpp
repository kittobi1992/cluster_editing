/*+++++++++++++++++ class EdgeList +++++++++++++++++++ */

#include <limits>

#include <vertexlists.h>

/* #################  Constructor  ##################### */


// empty constructor
VertexLists::VertexLists() : _size(0), _number_of_lists(0), nil(std::numeric_limits<unsigned short>::max()) {
}


// constructor which saves number of possible elements and number of possible lists
VertexLists::VertexLists(int size, pos_int_type number_of_lists) : _size(size), _number_of_lists(number_of_lists),
                                                                   nil(std::numeric_limits<unsigned short>::max()) {
    // create list of elements
    _vertex_list = vertex_list_type(_size);

    // save elements such that they belong to no lists
    for (int i = 0; i < _size; i++) {
        vertex_type v = {nil, nil, nil};
        _vertex_list[i] = v;
    };

    // create head list
    _heads = head_list_type(_number_of_lists);

    // save head elements as empty lists
    for (int i = 0; i < _number_of_lists; i++) {
        head_type hv = {nil, 0};
        _heads[i] = hv;
    };
}


/* #################  Public methods  ##################### */


// returns number of possible elements
int VertexLists::getSize() {
    return _size;
}


// returns number of lists
int VertexLists::getNumberOfLists() {
    return _number_of_lists;
}


// return list number of element i
VertexLists::pos_int_type VertexLists::getListNr(pos_int_type i) {
    return _vertex_list[i].list;
}


// return size of list "list"
int VertexLists::getListSize(pos_int_type list) {
    return _heads[list].size;
}


// return top element of list "list"
VertexLists::pos_int_type VertexLists::getTop(pos_int_type list) {
    if (_heads[list].size == 0) {
        throw VertexListsException("list is empty");
    }

    return _heads[list].head;
}


// return all elements in list "list" in a vertex_set_type list
VertexLists::vertex_set_type VertexLists::getList(pos_int_type list) {
    // get list size
    int size = _heads[list].size;

    // create vertex_set_type list
    vertex_set_type vertex_list = vertex_set_type(size);

    // set help pointer to first element
    pos_int_type help = _heads[list].head;

    // insert head element to output list
    vertex_list[0] = help;

    // iterate through internal list and add elements to output list
    for (int i = 1; i < size; i++) {
        help = _vertex_list[help].next;
        vertex_list[i] = help;
    }

    return vertex_list;
}


// add element i to list "list", just possible if i is not member of any list
void VertexLists::addToList(pos_int_type i, pos_int_type list) {
    // throw exception if i is already in other list
    if (getListNr(i) != nil) {
        throw VertexListsException("vertex already in another list");
    } else {
        // create element
        vertex_type v = {nil, nil, list};

        // if list not empty add it as first element
        if (_heads[list].size != 0) {
            pos_int_type oldRef = _heads[list].head;
            _vertex_list[oldRef].prev = i;
            v.next = oldRef;
        }

        // insert element in array
        _vertex_list[i] = v;

        // save inserted element as head
        _heads[list].head = i;
        _heads[list].size++;
    }
}


// deletes element i from its list, just possible if i is member of any list
void VertexLists::deleteFromList(pos_int_type i) {
    // get list number
    pos_int_type old_list = getListNr(i);

    // throw exception if i is in no list
    if (old_list == nil) {
        throw VertexListsException("vertex not in list");
    }

    // get element from array
    vertex_type v = _vertex_list[i];

    // set next element if existent to prev
    if (v.next != nil) {
        _vertex_list[v.next].prev = _vertex_list[i].prev;
    }

    // set pref element to next element
    if (i != _heads[old_list].head) {
        _vertex_list[v.prev].next = _vertex_list[i].next;
    } else {
        // set next element as head element
        if (v.next != nil) {
            _heads[old_list].head = v.next;
            // set list to empty
        } else {
            _heads[old_list].head = nil;
        }
    };

    // decrease list size
    _heads[old_list].size--;

    // set vertex i to blank
    vertex_type new_vertex = {nil, nil, nil};
    _vertex_list[i] = new_vertex;

}


// copy vertex i to list "list", just possible i is member of a list
void VertexLists::copyToList(pos_int_type i, pos_int_type list) {
    if (getListNr(i) != list) {
        deleteFromList(i);
        addToList(i, list);
    }
}


// merges two lists i,j, such that all elements from j are appended
// to list to i and list j is empty at its end
void VertexLists::mergeLists(pos_int_type i, pos_int_type j) {
    if (i != j && _heads[j].size != 0) {

        // set help pointer to head of j
        pos_int_type help = _heads[j].head;

        // iterate through elements of j
        // and set all elements "list" variable to list i
        while (_vertex_list[help].next != nil) {
            _vertex_list[help].list = i;
            help = _vertex_list[help].next;
        }
        _vertex_list[help].list = i;

        // set last element of j successor-variable to head of i
        _vertex_list[help].next = _heads[i].head;

        // if possible set original head of i precursor-variable to last element of j
        if (_heads[i].size != 0) {
            _vertex_list[_heads[i].head].prev = help;
        }

        // increase list i by elements of j and set new head of i
        _heads[i].size += _heads[j].size;
        _heads[i].head = _heads[j].head;

        // set list j to empty
        _heads[j].head = nil;
        _heads[j].size = 0;
    }

}


// print all lists on std output
void VertexLists::printLists() {
    for (int i = 0; i < _heads.size(); i++) {
        printList(i);
    }
}


// print list "list" on std output
void VertexLists::printList(pos_int_type list) {
    std::cout << "List " << list << " (size: " << _heads[list].size << " )";

    pos_int_type ref = _heads[list].head;
    while (ref != nil) {
        std::cout << "(" << ref << ") ";
        ref = _vertex_list[ref].next;
    }

    std::cout << std::endl;
}
