#include <limits>
#include <memory>
#include <algorithm>
#include <unordered_set>

#include "ShortestPathHeuristic.h"


ShortestPathHeuristic::ShortestPathHeuristic(size_t source, size_t graph_size, const AdjacencyMatrix &adj_matrix)
    : source(source), all_nodes(graph_size+1, nullptr) {
    size_t num_of_objectives = adj_matrix.get_num_of_objectives();
    size_t i = 0;
    for (auto node_iter = this->all_nodes.begin(); node_iter != this->all_nodes.end(); node_iter++) {
        *node_iter = std::make_shared<Node>(i++, std::vector<size_t>(num_of_objectives, 0), std::vector<size_t>(num_of_objectives, MAX_COST));
 }

    for (int j=0; j < num_of_objectives; j ++){
        compute(j, adj_matrix);
    }
}


std::vector<size_t> ShortestPathHeuristic::operator()(size_t node_id) {
    return this->all_nodes[node_id]->h;
}


// Implements Dijkstra shortest path algorithm per cost_idx cost function
void ShortestPathHeuristic::compute(size_t cost_idx, const AdjacencyMatrix &adj_matrix) {
    // Init all heuristics to MAX_COST
    for (auto node_iter = this->all_nodes.begin(); node_iter != this->all_nodes.end(); node_iter++) {
        (*node_iter)->h[cost_idx] = MAX_COST;
    }

    NodePtr node;
    NodePtr next;

    // Init open heap
    Node::more_than_specific_heurisitic_cost more_than(cost_idx);
    std::vector<NodePtr> open;
    std::make_heap(open.begin(), open.end(), more_than);

    this->all_nodes[this->source]->h[cost_idx] = 0;
    open.push_back(this->all_nodes[this->source]);
    std::push_heap(open.begin(), open.end(), more_than);


    while (open.empty() == false) {
        // Pop min from queue and process
        std::pop_heap(open.begin(), open.end(), more_than);
        node = open.back();
        open.pop_back();

        // Check to which neighbors we should extend the paths
        const std::vector<Edge> &outgoing_edges = adj_matrix[node->id];
        for(auto p_edge = outgoing_edges.begin(); p_edge != outgoing_edges.end(); p_edge++) {
            next = this->all_nodes[p_edge->target];

            // Dominance check
            if (next->h[cost_idx] <= (node->h[cost_idx]+p_edge->cost[cost_idx])) {
                continue;
            }

            // If not dominated push to queue
            next->h[cost_idx] = node->h[cost_idx] + p_edge->cost[cost_idx];
            open.push_back(next);
            std::push_heap(open.begin(), open.end(), more_than);
        }
    }
}


PartialShortestPathHeuristic::PartialShortestPathHeuristic(size_t source, size_t goal, const AdjacencyMatrix &adj_matrix)
    : source(source)  {
    num_obj = adj_matrix.get_num_of_objectives();

    for (int j=0; j < num_obj; j ++){
        data.push_back(std::unordered_map<size_t,size_t>());
        compute(j, adj_matrix, goal);
    }
    assert (default_values.size() == 2);
}




void PartialShortestPathHeuristic::set_value(size_t state, size_t cost_idx, size_t value){
    // if (data[cost_idx].find(state) == data.end()){
    //     data[cost_idx][state] = std::vector<size_t>(num_obj, 0);
    // }

    data[cost_idx][state] = value;
}


// Implements Dijkstra shortest path algorithm per cost_idx cost function
void PartialShortestPathHeuristic::compute(size_t cost_idx, const AdjacencyMatrix &adj_matrix, size_t goal) {
    // Init all heuristics to MAX_COST

    // Init open heap
    std::vector<std::pair<size_t, size_t>> open;
    // std::make_heap(open.begin(), open.end(), more_than);

    // this->all_nodes[this->source]->h[cost_idx] = 0;
    set_value(this->source, cost_idx, 0);
    open.push_back({0, this->source});
    auto order = std::greater<std::pair<size_t, size_t>>();
    std::push_heap(open.begin(), open.end(), order);

    std::unordered_map<size_t, size_t> generated;
    std::unordered_set<size_t> closed;
    size_t previous = 0;

    while (open.empty() == false) {
        // Pop min from queue and process
        std::pop_heap(open.begin(), open.end(), order);
        std::pair<size_t,size_t> curr = open.back();
        open.pop_back();
        size_t cost = curr.first;
        size_t node_id = curr.second;
        assert(cost >= previous);
        previous = cost;

        if (closed.find(node_id) != closed.end()){
            continue;
        }

        closed.insert(node_id);
        set_value(node_id, cost_idx, cost);

        if (node_id == goal){
            default_values.push_back(cost);
            break;
        }

        // Check to which neighbors we should extend the paths
        const std::vector<Edge> &outgoing_edges = adj_matrix[node_id];
        for(auto p_edge = outgoing_edges.begin(); p_edge != outgoing_edges.end(); p_edge++) {
            // Dominance check
            size_t ch_id = p_edge->target;
            if (generated.find(ch_id)!= generated.end() &&
                generated[ch_id]<= (cost+p_edge->apex[cost_idx])) {
                continue;
            }

            generated[ch_id]= (cost + p_edge->apex[cost_idx]);
            // If not dominated push to queue
            // next->h[cost_idx] = node->h[cost_idx] + p_edge->cost[cost_idx];
            open.push_back({cost + p_edge->apex[cost_idx], ch_id});
            std::push_heap(open.begin(), open.end(), order);
        }
    }
}




std::vector<size_t> PartialShortestPathHeuristic::operator()(size_t node_id) {
    auto v = default_values;
    for (size_t cost_idx = 0; cost_idx < v.size(); cost_idx++)
    if (data[cost_idx].find(node_id) != data[cost_idx].end()){
        v[cost_idx] = data[cost_idx][node_id];
    }
    return v;
}

