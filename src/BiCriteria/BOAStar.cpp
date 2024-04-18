#include <memory>
#include <algorithm>
#include <time.h>
#include <set>
#include <unordered_map>

#include "BOAStar.h"

BOAStar::BOAStar(const AdjacencyMatrix &adj_matrix, Pair<double> eps, const LoggerPtr logger) :
    AbstractSolver(adj_matrix, {eps[0], eps[1]}, logger) {}

void BOAStar::operator()(size_t source, size_t target, Heuristic &heuristic, SolutionSet &solutions, unsigned int time_limit) {
    // int time_limit = 300;
    start_time = std::clock();
    this->start_logging(source, target);

    NodePtr node;
    NodePtr next;

    // Saving all the unused NodePtrs in a vector improves performace for some reason
    std::vector<NodePtr> closed;

    // Vector to hold mininum cost of 2nd criteria per node
    std::vector<size_t> min_g2(this->adj_matrix.size()+1, MAX_COST);

    // Init open heap
    Node::more_than_full_cost more_than;
    std::vector<NodePtr> open;
    std::make_heap(open.begin(), open.end(), more_than);

    node = std::make_shared<Node>(source, std::vector<size_t>(2,0), heuristic(source));
    open.push_back(node);
    std::push_heap(open.begin(), open.end(), more_than);

    while (open.empty() == false) {
        if ((std::clock() - start_time)/CLOCKS_PER_SEC > time_limit){

            this->end_logging(solutions, false);
            return;
        }

        // Pop min from queue and process
        std::pop_heap(open.begin(), open.end(), more_than);
        node = open.back();
        open.pop_back();
        num_generation +=1;

        // Dominance check
        if ((((1+this->eps[1])*node->f[1]) >= min_g2[target]) ||
            (node->g[1] >= min_g2[node->id])) {
            closed.push_back(node);
            continue;
        }

        min_g2[node->id] = node->g[1];
        num_expansion += 1;


        if (node->id == target) {
            solutions.push_back(node);
            log_solution(node);
            continue;
        }

        // Check to which neighbors we should extend the paths
        const std::vector<Edge> &outgoing_edges = adj_matrix[node->id];
        for(auto p_edge = outgoing_edges.begin(); p_edge != outgoing_edges.end(); p_edge++) {
            size_t next_id = p_edge->target;
            std::vector<size_t> next_g = {node->g[0]+p_edge->cost[0], node->g[1]+p_edge->cost[1]};
            auto next_h = heuristic(next_id);

            // Dominance check
            if ((((1+this->eps[1])*(next_g[1]+next_h[1])) >= min_g2[target]) ||
                (next_g[1] >= min_g2[next_id])) {
                continue;
            }

            // If not dominated create node and push to queue
            // Creation is defered after dominance check as it is
            // relatively computational heavy and should be avoided if possible
            next = std::make_shared<Node>(next_id, next_g, next_h, node);

            open.push_back(next);
            std::push_heap(open.begin(), open.end(), more_than);

            closed.push_back(node);
        }
    }

    this->end_logging(solutions);
}


bool is_dominated(std::set<std::pair<int, int>>& pair_set, std::pair<int, int> vec){
    if ( pair_set.size() == 0 ){return false;}
    auto it = pair_set.lower_bound(vec);
    if (it != pair_set.end() && it->first == vec.first && it->second == vec.second){
        return true;
    }
    if (it != pair_set.begin()){
        auto it_prev = std::prev(it);
        if (it_prev->second <= vec.second){
            return true;
        }
    }
    return false;
}

void NAMOA::operator()(size_t source, size_t target, Heuristic &heuristic, SolutionSet &solutions, unsigned int time_limit) {
    // int time_limit = 300;
    start_time = std::clock();
    this->start_logging(source, target);

    NodePtr node;
    NodePtr next;

    // Saving all the unused NodePtrs in a vector improves performace for some reason
    std::vector<NodePtr> closed;

    // Vector to hold mininum cost of 2nd criteria per node
    // std::vector<size_t> min_g2(this->adj_matrix.size()+1, MAX_COST);
    std::unordered_map<int, std::set<std::pair<int, int>>> expanded_costs;

    // Init open heap
    Node::more_than_full_cost more_than;
    std::vector<NodePtr> open;
    std::make_heap(open.begin(), open.end(), more_than);
    int min_f2 = INT_MAX;

    node = std::make_shared<Node>(source, std::vector<size_t>(2,0), heuristic(source));
    open.push_back(node);
    std::push_heap(open.begin(), open.end(), more_than);

    while (open.empty() == false) {
        if ((std::clock() - start_time)/CLOCKS_PER_SEC > time_limit){

            this->end_logging(solutions, false);
            return;
        }

        // Pop min from queue and process
        std::pop_heap(open.begin(), open.end(), more_than);
        node = open.back();
        open.pop_back();
        num_generation +=1;

        // Dominance check
        if (((1+this->eps[1])*node->f[1]) >= min_f2){
            closed.push_back(node);
            continue;
        }

        if (node->id == target) {
            solutions.push_back(node);
            min_f2 = node->g[1];
            log_solution(node);
            continue;
        }

        if (expanded_costs.find(node->id) != expanded_costs.end() &&
            is_dominated(expanded_costs[node->id], {node->g[0], node->g[1]})) {
            closed.push_back(node);
            continue;
        }

        expanded_costs[node->id].insert({node->g[0], node->g[1]});
        num_expansion += 1;


        // Check to which neighbors we should extend the paths
        const std::vector<Edge> &outgoing_edges = adj_matrix[node->id];
        for(auto p_edge = outgoing_edges.begin(); p_edge != outgoing_edges.end(); p_edge++) {
            size_t next_id = p_edge->target;
            std::vector<size_t> next_g = {node->g[0]+p_edge->cost[0], node->g[1]+p_edge->cost[1]};
            auto next_h = heuristic(next_id);

            // Dominance check
            if ((((1+this->eps[1])*(next_g[1]+next_h[1])) >= min_f2) ||
                // (next_g[1] >= min_g2[next_id]
                is_dominated(expanded_costs[next_id], {next_g[0], next_g[1]})
                 ) {
                continue;
            }

            // If not dominated create node and push to queue
            // Creation is defered after dominance check as it is
            // relatively computational heavy and should be avoided if possible
            next = std::make_shared<Node>(next_id, next_g, next_h, node);

            open.push_back(next);
            std::push_heap(open.begin(), open.end(), more_than);

            closed.push_back(node);
        }
    }

    this->end_logging(solutions);
}



bool is_dom(const std::pair<int, int>& p1, const std::pair<int, int>& p2){
    return p1.first <= p2.first && p1.second <= p2.second;
};

void remove_dominated(std::set<std::pair<int,int>> & pair_set, std::pair<int,int> p){
    auto it = pair_set.upper_bound(p);
    if (it != pair_set.end()){
        auto curr = it;
        while (curr != pair_set.end() && is_dom(p, *curr)){
            curr = next(curr);
        }
        pair_set.erase(it, curr);
    }
}

void BOAStarMVH::operator()(size_t source, size_t target, Heuristic &heuristic, SolutionSet &solutions, unsigned int time_limit) {
    // int time_limit = 300;
    start_time = std::clock();
    this->start_logging(source, target);

    NodePtr node;
    NodePtr next;

    // Saving all the unused NodePtrs in a vector improves performace for some reason
    std::vector<NodePtr> closed;

    // Vector to hold mininum cost of 2nd criteria per node
    // std::vector<size_t> min_g2(this->adj_matrix.size()+1, MAX_COST);
    std::unordered_map<int, std::set<std::pair<int, int>>> expanded_costs;

    // Init open heap
    Node::more_than_full_cost more_than;
    std::vector<NodePtr> open;
    std::make_heap(open.begin(), open.end(), more_than);
    int min_f2 = INT_MAX;

    node = std::make_shared<Node>(source, std::vector<size_t>(2,0), std::vector<size_t>(2,0));
    open.push_back(node);
    std::push_heap(open.begin(), open.end(), more_than);

    while (open.empty() == false) {
        if ((std::clock() - start_time)/CLOCKS_PER_SEC > time_limit){

            this->end_logging(solutions, false);
            return;
        }

        // Pop min from queue and process
        std::pop_heap(open.begin(), open.end(), more_than);
        node = open.back();
        open.pop_back();
        num_generation +=1;

        // Dominance check
        if (((1+this->eps[1])*node->f[1]) >= min_f2){
            auto new_h = mvh->heuristic(node->id, node->g, min_f2);
            if (new_h.second + node->g[1] < min_f2){
                node->update_h({ (size_t) new_h.first, (size_t) new_h.second});
                open.push_back(node);
                std::push_heap(open.begin(), open.end(), more_than);
            } else {
                closed.push_back(node);
            }
            continue;
        }

        if (node->id == target) {
            solutions.push_back(node);
            min_f2 = node->g[1];
            log_solution(node);
            continue;
        }

        if (expanded_costs.find(node->id) != expanded_costs.end() &&
            is_dominated(expanded_costs[node->id], {node->g[0], node->g[1]})) {
            closed.push_back(node);
            continue;
        }

        expanded_costs[node->id].insert({node->g[0], node->g[1]});
        num_expansion += 1;


        // Check to which neighbors we should extend the paths
        const std::vector<Edge> &outgoing_edges = adj_matrix[node->id];
        for(auto p_edge = outgoing_edges.begin(); p_edge != outgoing_edges.end(); p_edge++) {
            size_t next_id = p_edge->target;
            std::vector<size_t> next_g = {node->g[0]+p_edge->cost[0], node->g[1]+p_edge->cost[1]};
            // auto next_h = heuristic(next_id);
            auto h_p = mvh->heuristic(next_id, next_g, min_f2);
            std::vector<size_t> next_h = {h_p.first, h_p.second};

            // Dominance check
            if ((((1+this->eps[1])*(next_g[1]+next_h[1])) >= min_f2) ||
                // (next_g[1] >= min_g2[next_id]
                is_dominated(expanded_costs[next_id], {next_g[0], next_g[1]})
                 ) {
                continue;
            }

            // If not dominated create node and push to queue
            // Creation is defered after dominance check as it is
            // relatively computational heavy and should be avoided if possible
            next = std::make_shared<Node>(next_id, next_g, next_h, node);

            open.push_back(next);
            std::push_heap(open.begin(), open.end(), more_than);

            closed.push_back(node);
        }
    }

    this->end_logging(solutions);
}

inline bool is_dominated(SolutionSet & solutions, NodePtr node, Pair<double> eps = {0,0}){
    for (auto sol: solutions){
        if (sol->g[0] <= (1 + eps[0]) * node->g[0] + node->h[0] &&
            sol->g[1] <= (1 + eps[1]) * node->g[1] + node->h[1]
            ){
            return true;
        }
    }
    return false;
}




void BOAStar::log_solution(NodePtr node){
    solution_log.push_back({std::clock() - start_time, node});
}
