#include <algorithm>
#include "Utils/MapQueue.h"

template<class Node, class more_than_full_cost>
MapQueueHash<Node, more_than_full_cost>::MapQueueHash(){

    std::make_heap(this->heap.begin(), this->heap.end(), this->more_than);
}

template<class Node, class more_than_full_cost>
bool MapQueueHash<Node, more_than_full_cost>::empty() {
    return this->heap.empty();
}


template<class Node, class more_than_full_cost>
Node MapQueueHash<Node, more_than_full_cost>::pop() {
    // Pop from min heap
    std::pop_heap(this->heap.begin(), this->heap.end(), this->more_than);
    Node pp = this->heap.back();
    this->heap.pop_back();

    // Remove from open map
    std::list<Node> &relevant_pps = this->open_map[pp->id];
    for (auto iter = relevant_pps.begin(); iter != relevant_pps.end(); ++iter) {
        if (pp == *iter) {
            relevant_pps.erase(iter);
            break;
        }
    }

    return pp;
}

template<class Node, class more_than_full_cost>
void MapQueueHash<Node, more_than_full_cost>::insert(Node &pp) {
    // Insert to min heap
    this->heap.push_back(pp);
    std::push_heap(this->heap.begin(), this->heap.end(), this->more_than);

    // Insert to open map
    this->open_map[pp->id].push_back(pp);
}

template<class Node, class more_than_full_cost>
std::list<Node> &MapQueueHash<Node, more_than_full_cost>::get_open(size_t id) {
	return this->open_map[id];
}


template class MapQueueHash<ApexPathPairPtr, ApexPathPair::more_than_full_cost>;
