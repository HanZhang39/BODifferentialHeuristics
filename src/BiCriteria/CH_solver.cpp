#include "CH_solver.h"
#include "Utils/IOUtils.h"
#include <fstream>


CHGraphRT::CHGraphRT(std::string fname, size_t num_obj){
    this->num_obj = num_obj;

    std::ifstream  file(fname.c_str());

    // assert (file.is_open() == false);


    std::string line;

    size_t graph_size, contraction_order_size;
    std::getline(file, line);

    std::vector<std::string> decomposed_line_1;
    split_string(line, "\t", decomposed_line_1);
    if (decomposed_line_1.size() == 1){
        std::cout << "loading old version file...\n";
        graph_size = std::stoul(decomposed_line_1[0]);
        contraction_order_size = std::stoul(decomposed_line_1[0]);
    } else {
        graph_size = std::stoul(decomposed_line_1[1]);
        contraction_order_size = std::stoul(decomposed_line_1[0]);
    }


    std::cout << "graph size:" << graph_size << std::endl;
    std::cout << "contract size:" << contraction_order_size << std::endl;

    state_to_order.resize(graph_size, graph_size + 10);

    for (size_t i = 0; i < contraction_order_size; i++){
        std::getline(file, line);

        if (line == "") {
            break;
        } else if (line[0] == '#') {
            continue; // Commented out queries
        }

        std::vector<std::string> decomposed_line;
        split_string(line, "\t", decomposed_line);
        size_t node_id = std::stoul(decomposed_line[0]);
        // if (node_id >=state_to_order.size()){
        //     std::cout << "error1" << std::endl;
        // }
        // if (state_to_order[node_id] != 0){
        //     std::cout << "errors" << std::endl;
        // }
        state_to_order[node_id] = i + 1;

    }

    std::cout << state_to_order.size() << std::endl;

    std::getline(file, line);
    size_t shortcut_num = std::stoul(line);
    std::cout << "shortcuts:" << shortcut_num << std::endl;


    std::vector<Edge> up_edges;
    std::vector<Edge> down_edges;
    std::vector<Edge> equal_edges;

    for (size_t i = 0; i < shortcut_num; i++){
        std::getline(file, line);

        if (line == "") {
            break;
        } else if (line[0] == '#') {
            continue; // Commented out queries
        }

        std::vector<std::string> decomposed_line;
        split_string(line, "\t", decomposed_line);
        if (decomposed_line.size() < 6){
            std::cout << "errors" << i << std::endl;
        }
        size_t state_from = std::stoul(decomposed_line[0]);
        size_t state_to = std::stoul(decomposed_line[1]);

        std::vector<size_t> cost;
        std::vector<size_t> apex;
        for (size_t j = 0 ; j < num_obj; j++){
            cost.push_back(std::stoul(decomposed_line[2 + j]));
        }
        for (size_t j = 0 ; j < num_obj; j++){
            apex.push_back(std::stoul(decomposed_line[2 + num_obj + j]));
        }
        Edge e(state_from, state_to, cost, apex);
        if (state_to_order[state_to] > state_to_order[state_from]){
            up_edges.push_back(e);
        } else if (state_to_order[state_to] < state_to_order[state_from]){
            down_edges.push_back(e);
        } else{
            equal_edges.push_back(e);
            // std::cout << "error!" << std::endl;
        }
    }
    up_graph = AdjacencyMatrix(graph_size, up_edges);
    up_inv_graph = AdjacencyMatrix(graph_size, up_edges, true);

    down_graph = AdjacencyMatrix(graph_size, down_edges);
    down_inv_graph = AdjacencyMatrix(graph_size, down_edges, true);

    if (equal_edges.size() != 0){
        std::cout << "equal weight edge: " << equal_edges.size() << std::endl;
    }
    for (const auto& edge:equal_edges){
        up_graph.add(edge);
        down_inv_graph.add(edge.inverse());
        // down_graph.add(edge);
        // up_inv_graph.add(edge.inverse());
    }

}



std::unordered_set<size_t> CHGraphRT::all_up_states(size_t state){
    std::unordered_set<size_t> closed;

    std::list<size_t> open;
    open.push_back(state);
    closed.insert(state);

    while (!open.empty()){
        size_t curr = open.back();
        open.pop_back();

        for (const auto&e:up_graph[curr]){
            if ( closed.find(e.target) == closed.end() ){
                closed.insert(e.target);
                open.push_back(e.target);
            }
        }
    }

    return closed;
}
std::unordered_set<size_t> CHGraphRT::all_up_states_r(size_t state){
    std::unordered_set<size_t> closed;

    std::list<size_t> open;
    open.push_back(state);
    closed.insert(state);

    while (!open.empty()){
        size_t curr = open.back();
        open.pop_back();

        for (const auto&e : down_inv_graph[curr]){
            if ( closed.find(e.target) == closed.end() ){
                closed.insert(e.target);
                open.push_back(e.target);
            }
        }
    }
    return closed;
}


CHShortestPathHeuristic::CHShortestPathHeuristic(size_t source,
                                                 const CHGraphRT& chg,
                                                 const std::unordered_set<size_t>& up_states,
                                                 const std::unordered_set<size_t>& down_states)
    : source(source)  {
    num_obj = chg.get_num_of_objectives();

    for (size_t j=0; j < num_obj; j ++){
        data.push_back(std::unordered_map<size_t,size_t>());
        compute(j, chg, up_states, down_states);
    }
}



// Implements Dijkstra shortest path algorithm per cost_idx cost function
void CHShortestPathHeuristic::compute(size_t cost_idx, const CHGraphRT& chg, const std::unordered_set<size_t>& up_states, const std::unordered_set<size_t>& down_states){
    // Init all heuristics to MAX_COST

    // Init open heap
    std::vector<std::pair<size_t, size_t>> open;
    data[cost_idx][source] = 0;
    open.push_back({0, source});
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
        data[cost_idx][node_id] = cost;

        // Check to which neighbors we should extend the paths
        if (down_states.find(node_id) != down_states.end()){
            const std::vector<Edge> &outgoing_edges = chg.down_inv_edges(node_id);
            for(auto p_edge = outgoing_edges.begin(); p_edge != outgoing_edges.end(); p_edge++) {
                // Dominance check
                size_t ch_id = p_edge->target;
                if (up_states.find(ch_id) == up_states.end() && down_states.find(ch_id) == down_states.end()){
                    continue;
                }
                if (generated.find(ch_id)!= generated.end() &&
                    generated[ch_id]<= (cost+p_edge->apex[cost_idx])) {
                    continue;
                }

                generated[ch_id]= (cost + p_edge->apex[cost_idx]);
                open.push_back({cost + p_edge->apex[cost_idx], ch_id});
                std::push_heap(open.begin(), open.end(), order);
            }
        }
        if (up_states.find(node_id) != up_states.end()){
            const std::vector<Edge> &outgoing_edges = chg.up_inv_edges(node_id);
            for(auto p_edge = outgoing_edges.begin(); p_edge != outgoing_edges.end(); p_edge++) {
                // Dominance check
                size_t ch_id = p_edge->target;
                if (up_states.find(ch_id) == up_states.end() && down_states.find(ch_id) == down_states.end()){
                    continue;
                }
                if (generated.find(ch_id)!= generated.end() &&
                    generated[ch_id]<= (cost+p_edge->apex[cost_idx])) {
                    continue;
                }

                generated[ch_id]= (cost + p_edge->apex[cost_idx]);
                open.push_back({cost + p_edge->apex[cost_idx], ch_id});
                std::push_heap(open.begin(), open.end(), order);
            }
        }
    }
}

std::vector<size_t> CHShortestPathHeuristic::operator()(size_t node_id) {
    auto v = std::vector<size_t>(data.size(), 0);
    for (size_t cost_idx = 0; cost_idx < v.size(); cost_idx++)
    if (data[cost_idx].find(node_id) != data[cost_idx].end()){
        v[cost_idx] = data[cost_idx][node_id];
    }
    return v;
}


void BOAStarCH::operator()(size_t source, size_t target, SolutionSet &solutions, unsigned int time_limit) {
    start_time = std::clock();
    auto up_states = chg.all_up_states(source);
    auto down_states = chg.all_up_states_r(target);
    if (source != 0 and target != 0){
        std::cout << up_states.size() << std::endl;
        std::cout << down_states.size() << std::endl;
    }

    auto heuristic = CHShortestPathHeuristic(target, chg, up_states, down_states);

    for (auto &l: heuristic.data){
        std::cout << l.size() << std::endl;
    }

    preprocess_time = std::clock() - start_time;

    NodePtr node;
    NodePtr next;

    // Saving all the unused NodePtrs in a vector improves performace for some reason
    std::vector<NodePtr> closed;

    // Vector to hold mininum cost of 2nd criteria per node
    std::unordered_map<size_t, size_t> min_g2;

    // Init open heap
    Node::more_than_full_cost more_than;
    std::vector<NodePtr> open;
    std::make_heap(open.begin(), open.end(), more_than);

    node = std::make_shared<Node>(source, std::vector<size_t>(2,0), heuristic(source));
    open.push_back(node);
    std::push_heap(open.begin(), open.end(), more_than);

    std::cout <<"hinit "<< heuristic(source)[0] << ", "<< heuristic(source)[1] << std::endl;

    while (open.empty() == false) {
        if ((std::clock() - start_time)/CLOCKS_PER_SEC > time_limit){

            total_time = std::clock() - start_time;
            return;
        }

        // Pop min from queue and process
        std::pop_heap(open.begin(), open.end(), more_than);
        node = open.back();
        open.pop_back();
        num_generation +=1;

        // Dominance check
        if (min_g2.find(target) != min_g2.end()  &&
            node->f[1] >= min_g2[target]){
            continue;
        }
        if (min_g2.find(node->id) != min_g2.end()  &&
            node->g[1] >= min_g2[node->id]){
            continue;
        }

        min_g2[node->id] = node->g[1];
        num_expansion += 1;


        if (node->id == target) {
            solutions.push_back(node);
            continue;
        }

        // Check to which neighbors we should extend the paths
        if (up_states.find(node->id) != up_states.end()){
            const std::vector<Edge> &outgoing_edges = chg.up_edges(node->id);
        for(auto p_edge = outgoing_edges.begin(); p_edge != outgoing_edges.end(); p_edge++) {
            size_t next_id = p_edge->target;
            if (up_states.find(next_id) == up_states.end() && down_states.find(next_id) == down_states.end()){
                continue;
            }
            std::vector<size_t> next_g = {node->g[0]+p_edge->cost[0], node->g[1]+p_edge->cost[1]};
            auto next_h = heuristic(next_id);

            // Dominance check
            if (min_g2.find(target) != min_g2.end()  &&
                next_g[1] + next_h[1] >= min_g2[target]){
                continue;
            }
            if (min_g2.find(next_id) != min_g2.end()  &&
                next_g[1] >= min_g2[next_id]){
                continue;
            }

            next = std::make_shared<Node>(next_id, next_g, next_h, node);

            open.push_back(next);
            std::push_heap(open.begin(), open.end(), more_than);
        }
        }
        if (down_states.find(node->id) != down_states.end()){
            const std::vector<Edge> &outgoing_edges = chg.down_edges(node->id);
            for(auto p_edge = outgoing_edges.begin(); p_edge != outgoing_edges.end(); p_edge++) {
                size_t next_id = p_edge->target;
                if (up_states.find(next_id) == up_states.end() && down_states.find(next_id) == down_states.end()){
                    continue;
                }
                std::vector<size_t> next_g = {node->g[0]+p_edge->cost[0], node->g[1]+p_edge->cost[1]};
                auto next_h = heuristic(next_id);

                // Dominance check
                if (min_g2.find(target) != min_g2.end()  &&
                    next_g[1] + next_h[1] >= min_g2[target]){
                    continue;
                }
                if (min_g2.find(next_id) != min_g2.end()  &&
                    next_g[1] >= min_g2[next_id]){
                    continue;
                }

                next = std::make_shared<Node>(next_id, next_g, next_h, node);

                open.push_back(next);
                std::push_heap(open.begin(), open.end(), more_than);
            }
        }
    }

    total_time = std::clock() - start_time;


}



void ApexSearchCH::operator()(size_t source, size_t target, SolutionSet &solutions, unsigned int time_limit) {
    start_time = std::clock();
    auto up_states = chg.all_up_states(source);
    auto down_states = chg.all_up_states_r(target);
    if (source != 0 and target != 0){
        std::cout << up_states.size() << std::endl;
        std::cout << down_states.size() << std::endl;
    }

    auto sp_heuristic = CHShortestPathHeuristic(target, chg, up_states, down_states);

    using std::placeholders::_1;
    Heuristic heuristic = std::bind( &CHShortestPathHeuristic::operator(), sp_heuristic, _1);

    for (auto &l: sp_heuristic.data){
        std::cout << l.size() << std::endl;
    }

    preprocess_time = std::clock() - start_time;
    size_t num_of_objectives = chg.get_num_of_objectives();
    if (num_of_objectives == 2){
        solution_dom_checker = std::make_unique<SolutionCheck>(eps);
    }else{
        solution_dom_checker = std::make_unique<SolutionCheckLinear>(eps);
    }

    ApexPathSolutionSet ap_solutions;
    ApexPathPairPtr   ap;
    ApexPathPairPtr   next_ap;

    // Saving all the unused PathPairPtrs in a vector improves performace for some reason
    // std::vector<ApexPathPairPtr> closed;

    // Vector to hold mininum cost of 2nd criteria per node
    // std::vector<size_t> min_g2(this->adj_matrix.size()+1, MAX_COST);
    
    // Init open heap
    APQueueHash open;

    NodePtr source_node = std::make_shared<Node>(source, std::vector<size_t>(num_of_objectives, 0), heuristic(source));
    ap = std::make_shared<ApexPathPair>(source_node, source_node, heuristic);
    open.insert(ap);

    while (open.empty() == false) {
        if ((std::clock() - start_time)/CLOCKS_PER_SEC > time_limit){

            return;
        }

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

        min_g2[ap->id] = ap->apex->g[1];

        num_expansion += 1;

        // expanded[ap->id].push_back(ap);

        if (ap->id == target) {
            this->merge_to_solutions(ap, ap_solutions);
            continue;
        }


        // Check to which neighbors we should extend the paths
        if (up_states.find(ap->id) != up_states.end()){
            const std::vector<Edge> &outgoing_edges = chg.up_edges(ap->id);
            for(auto p_edge = outgoing_edges.begin(); p_edge != outgoing_edges.end(); p_edge++) {
                size_t next_id = p_edge->target;
                if (up_states.find(next_id) == up_states.end() && down_states.find(next_id) == down_states.end()){
                    continue;
                }
                next_ap = std::make_shared<ApexPathPair>(ap, *p_edge);
                if (is_dominated(next_ap)){
                    continue;
                }
                this->insert(next_ap, open);
            }
        }
        if (down_states.find(ap->id) != down_states.end()){
            const std::vector<Edge> &outgoing_edges = chg.down_edges(ap->id);
            for(auto p_edge = outgoing_edges.begin(); p_edge != outgoing_edges.end(); p_edge++) {
                size_t next_id = p_edge->target;
                if (up_states.find(next_id) == up_states.end() && down_states.find(next_id) == down_states.end()){
                    continue;
                }
                next_ap = std::make_shared<ApexPathPair>(ap, *p_edge);
                if (is_dominated(next_ap)){
                    continue;
                }
                this->insert(next_ap, open);
            }
        }
    }

    total_time = std::clock() - start_time;

    for (auto solution = ap_solutions.begin(); solution != ap_solutions.end(); ++solution) {
        solutions.push_back((*solution)->path_node);

    }

}


void ApexSearchCH::insert(ApexPathPairPtr &ap, APQueueHash &queue){
    
    std::list<ApexPathPairPtr> &relevant_aps = queue.get_open(ap->id);
    for (auto existing_ap = relevant_aps.begin(); existing_ap != relevant_aps.end(); ++existing_ap) {
        if ((*existing_ap)->is_active == false) {
            continue;
        }
        if (ap->update_nodes_by_merge_if_bounded_g(*existing_ap, this->eps, ms) == true) {
            if ((ap-> apex!= (*existing_ap)->apex) ||
                (ap-> path_node!= (*existing_ap)->path_node)) {
                (*existing_ap)->is_active = false;
                queue.insert(ap);
            }
            return;
        }
    }
    queue.insert(ap);
}

bool ApexSearchCH::is_dominated(ApexPathPairPtr ap){
    if (min_g2.find(ap->id) != min_g2.end() && ap->apex->g[1] >= min_g2[ap->id]){
        return true;
    }


    return solution_dom_checker->is_dominated(ap);
}

void ApexSearchCH::merge_to_solutions(const ApexPathPairPtr &ap, ApexPathSolutionSet &solutions){
    for (auto existing_solution = solutions.begin(); existing_solution != solutions.end(); ++existing_solution) {
        if ((*existing_solution)->update_nodes_by_merge_if_bounded_g(ap, this->eps, ms) == true) {
            return;
        }
    }
    solutions.push_back(ap);
    // std::cout << "update solution checker" << std::endl;
    solution_dom_checker->add_node(ap);
}
