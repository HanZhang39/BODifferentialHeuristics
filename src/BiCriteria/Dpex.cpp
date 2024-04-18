#include <memory>
#include <vector>

#include <iostream>

#include "ApexSearch.h"



bool Dpex::is_dominated(ApexPathPairPtr ap){
    if (local_dom_checker->is_dominated(ap)){
        return true;
    }
    return false;
}


void Dpex::operator()(size_t source, std::vector<std::vector<std::pair<size_t, size_t>>> &solutions, unsigned int time_limit) {

    ZeroH zh(num_of_objectives);

    using std::placeholders::_1;
    Heuristic heuristic = std::bind( &ZeroH::get_val, zh, _1);

    solutions.clear();
    solutions.resize(adj_matrix.size() + 1);

    init_search();

    auto start_time = std::clock();

    if (num_of_objectives == 2){
        local_dom_checker = std::make_unique<LocalCheck>(eps, this->adj_matrix.size());
    }else{
        local_dom_checker = std::make_unique<LocalCheckLinear>(eps, this->adj_matrix.size());
    }


    // this->start_logging(source, target);

    ApexPathSolutionSet ap_solutions;
    ApexPathPairPtr   ap;
    ApexPathPairPtr   next_ap;

    // Saving all the unused PathPairPtrs in a vector improves performace for some reason
    // std::vector<ApexPathPairPtr> closed;

    // Vector to hold mininum cost of 2nd criteria per node
    // std::vector<size_t> min_g2(this->adj_matrix.size()+1, MAX_COST);
    
    // Init open heap
    APQueue open(this->adj_matrix.size()+1);

    NodePtr source_node = std::make_shared<Node>(source, std::vector<size_t>(num_of_objectives, 0),  std::vector<size_t>(num_of_objectives, 0));
    ap = std::make_shared<ApexPathPair>(source_node, source_node, heuristic);
    open.insert(ap);

    while (open.empty() == false) {
        if ((std::clock() - start_time)/CLOCKS_PER_SEC > time_limit){
            //            this->end_logging(solutions, false);
            return;
        }
        // Pop min from queue and process
        ap = open.pop();
        num_generation +=1;

        // Optimization: PathPairs are being deactivated instead of being removed so we skip them.
        if (ap->is_active == false) {
            continue;
        }

        // Dominance check
        if (is_dominated(ap)){
            continue;
        }

        //  min_g2[ap->id] = ap->bottom_right->g[1];
        local_dom_checker->add_node(ap);
        solutions[ap->id].push_back({ap->apex->g[0], ap->apex->g[1]});

        num_expansion += 1;

        expanded[ap->id].push_back(ap);

        // Check to which neighbors we should extend the paths
        const std::vector<Edge> &outgoing_edges = adj_matrix[ap->id];
        for(auto p_edge = outgoing_edges.begin(); p_edge != outgoing_edges.end(); p_edge++) {
            // Prepare extension of path pair

            next_ap = std::make_shared<ApexPathPair>(ap, *p_edge);

            // Dominance check
            // if ((((1+this->eps[1])*(bottom_right_next_g[1]+next_h[1])) >= min_g2[target]) ||
            //     (bottom_right_next_g[1] >= min_g2[next_id])) {
            if (is_dominated(next_ap)){
                continue;
            }

            // If not dominated extend path pair and push to queue
            // Creation is defered after dominance check as it is
            // relatively computational heavy and should be avoided if possible
            // std::cout <<"generate node on " << next_ap->id << std::endl;
            this->insert(next_ap, open);
            // closed.push_back(pp);
        }
    }

    // Pair solutions is used only for logging, as we need both the solutions for testing reasons

    //    this->end_logging(solutions);
}


void Dpex::merge_to_solutions(const ApexPathPairPtr &ap, ApexPathSolutionSet &solutions) {
    for (auto existing_solution = solutions.rbegin(); existing_solution != solutions.rend(); ++existing_solution) {
    if ((*existing_solution)->update_nodes_by_merge_if_bounded(ap, this->eps, ms) == true) {
        // std::cout << "merging" << std::endl;
        return;
    }
    }
    solutions.push_back(ap);
}

void Dpex::operator()(size_t source, std::unordered_map<size_t, std::vector<ApexPathPairPtr>> &solutions, unsigned int time_limit){
    
    ZeroH zh(num_of_objectives);

    using std::placeholders::_1;
    Heuristic heuristic = std::bind( &ZeroH::get_val, zh, _1);

    // init_search();
    AbstractSolver::init_search();


    auto start_time = std::clock();

    if (num_of_objectives == 2){
        local_dom_checker = std::make_unique<LocalCheck>(eps, this->adj_matrix.size());
    }else{
        local_dom_checker = std::make_unique<LocalCheckLinear>(eps, this->adj_matrix.size());
    }

    ApexPathPairPtr   ap;
    ApexPathPairPtr   next_ap;

    
    // Init open heap
    APQueue open(this->adj_matrix.size()+1);

    NodePtr source_node = std::make_shared<Node>(source, std::vector<size_t>(num_of_objectives, 0),  std::vector<size_t>(num_of_objectives, 0));
    ap = std::make_shared<ApexPathPair>(source_node, source_node, heuristic);
    open.insert(ap);

    while (open.empty() == false) {
        if ((std::clock() - start_time)/CLOCKS_PER_SEC > time_limit){
            //            this->end_logging(solutions, false);
            return;
        }
        // Pop min from queue and process
        ap = open.pop();
        num_generation +=1;

        // Optimization: PathPairs are being deactivated instead of being removed so we skip them.
        if (ap->is_active == false) {
            continue;
        }

        // Dominance check
        if (is_dominated(ap)){
            continue;
        }

        local_dom_checker->add_node(ap);
        if (solutions.find(ap->id) != solutions.end()){
            merge_to_solutions(ap, solutions[ap->id]);
        }

        num_expansion += 1;

        // Check to which neighbors we should extend the paths
        const std::vector<Edge> &outgoing_edges = adj_matrix[ap->id];
        for(auto p_edge = outgoing_edges.begin(); p_edge != outgoing_edges.end(); p_edge++) {
            // Prepare extension of path pair

            next_ap = std::make_shared<ApexPathPair>(ap, *p_edge);
            next_ap->parent = nullptr;

            if (is_dominated(next_ap)){
                continue;
            }

            this->insert(next_ap, open);
        }
    }

}
