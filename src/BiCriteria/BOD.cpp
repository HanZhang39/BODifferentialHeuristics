#include <memory>
#include <algorithm>
#include <time.h>

#include "BOAStar.h"


void BOD::operator()(size_t source, std::vector<std::vector<std::pair<size_t, size_t>>> &solutions, unsigned int time_limit) {
    ZeroH zh(2);
    using std::placeholders::_1;
    Heuristic heuristic = std::bind( &ZeroH::get_val, zh, _1);

    // int time_limit = 300;
    start_time = std::clock();

    solutions.clear();
    solutions.resize(adj_matrix.size() + 1);

    NodePtr node;
    NodePtr next;

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

            return;
        }

        // Pop min from queue and process
        std::pop_heap(open.begin(), open.end(), more_than);
        node = open.back();
        open.pop_back();
        num_generation +=1;

        // Dominance check
        if ((node->g[1] >= min_g2[node->id])) {
            continue;
        }

        min_g2[node->id] = node->g[1];
        solutions[node->id].push_back({ (size_t) node->g[0], (size_t) node->g[1] });
        num_expansion += 1;


        // Check to which neighbors we should extend the paths
        const std::vector<Edge> &outgoing_edges = adj_matrix[node->id];
        for(auto p_edge = outgoing_edges.begin(); p_edge != outgoing_edges.end(); p_edge++) {
            size_t next_id = p_edge->target;
            std::vector<size_t> next_g = {node->g[0]+p_edge->cost[0], node->g[1]+p_edge->cost[1]};
            auto next_h = heuristic(next_id);

            // Dominance check
            if(next_g[1] >= min_g2[next_id]) {
                continue;
            }

            // If not dominated create node and push to queue
            // Creation is defered after dominance check as it is
            // relatively computational heavy and should be avoided if possible
            next = std::make_shared<Node>(next_id, next_g, next_h, nullptr);

            open.push_back(next);
            std::push_heap(open.begin(), open.end(), more_than);

        }
    }

}
