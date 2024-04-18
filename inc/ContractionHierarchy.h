#include "Utils/Definitions.h"
#include "NAMOA.h"

class CHState{
public:
    size_t state;

    size_t delet_edge_num = 0;
    size_t add_edge_num = 0;
    size_t height = 0;
    double priority = 0;

    std::vector<Edge> shortcuts;
    bool flag=false;

    CHState(size_t state): state(state) {};

    struct smaller_weight {
        bool operator()(const CHState &a, const CHState &b) const;
    };

};

using CHStatePtr   = std::shared_ptr<CHState>;


struct ContractionOrdering{
    bool operator()(const CHStatePtr &a, const CHStatePtr &b) const;
};


class ApexContractionHierarchy{
private:
    clock_t total_time;
    std::vector<size_t> contraction_order;

    void update_priority(size_t state);

    size_t verbal=2;

    std::vector<double> epsilon;

    AdjacencyMatrix graph;
    AdjacencyMatrix inv_graph;

    std::vector<Edge> all_shortcuts;
    std::vector<CHStatePtr> states;

    // add shortcuts
    // void add_edge();

    void delete_state(size_t state);
    void contract_state(size_t state);
    void witness_search(size_t source, size_t target, std::vector<Edge> &edges) ;
    void witness_search3(size_t source, size_t target, std::vector<Edge> &edges) ;

    // update the state, unset its flag etc.
    void edge_difference(size_t state);
    void merge_edges(std::vector<Edge> & edges);



    // if found a edge that can be merged
    // return True, update new_edge, and update the ptr
    bool can_merge(Edge& new_edge, int & index);


public:

    clock_t get_time() {return total_time;}
    clock_t get_num_shortcuts() {return all_shortcuts.size();}

    ApexContractionHierarchy(std::vector<double> epsilon, size_t graph_size, std::vector<Edge> & edges): epsilon(epsilon), graph(graph_size, edges), inv_graph(graph_size, edges, true){};

    void contract(size_t contract_limit=0);

    void write_to(std::string fname);

};
