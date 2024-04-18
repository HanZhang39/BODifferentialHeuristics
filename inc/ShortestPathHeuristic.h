#ifndef EXAMPLE_SHORTEST_PATH_HEURISTIC_H
#define EXAMPLE_SHORTEST_PATH_HEURISTIC_H

#include "Utils/Definitions.h"
#include <unordered_map> 

// Precalculates heuristic based on Dijkstra shortest paths algorithm.
// On call to operator() returns the value of the heuristic in O(1)
class ShortestPathHeuristic {
private:
    size_t                  source;
    std::vector<NodePtr>    all_nodes;

    void compute(size_t cost_idx, const AdjacencyMatrix& adj_matrix);
public:
    ShortestPathHeuristic(size_t source, size_t graph_size, const AdjacencyMatrix &adj_matrix);
    std::vector<size_t> operator()(size_t node_id);
    void set_all_to_zero(){
        for (auto n: all_nodes){
            n->h = {0, 0};
        }
    }
};



class PartialShortestPathHeuristic {
private:
    size_t                  source;
    size_t num_obj;
    std::vector<size_t> default_values;

    std::vector<std::unordered_map<size_t, size_t>> data;

    void set_value(size_t state, size_t cost_idx, size_t value);
    void compute(size_t cost_idx, const AdjacencyMatrix& adj_matrix, size_t goal);
public:
    PartialShortestPathHeuristic(size_t source, size_t goal, const AdjacencyMatrix &adj_matrix);
    std::vector<size_t> operator()(size_t node_id);
};



#endif // EXAMPLE_SHORTEST_PATH_HEURISTIC_H
