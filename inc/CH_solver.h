#include "Utils/Definitions.h"
#include "ContractionHierarchy.h"
#include <unordered_set>
#include <unordered_map>
#include "DominanceChecker.h"
#include "Utils/MapQueue.h"

class CHGraphRT{
    // CH graph for runtime
    size_t num_obj= 0;

    AdjacencyMatrix up_graph;
    AdjacencyMatrix up_inv_graph;

    AdjacencyMatrix down_graph;
    AdjacencyMatrix down_inv_graph;

    std::vector<Edge> all_edges;


    // std::vector<size_t> contraction_order;
    std::vector<size_t> state_to_order;

public:

    std::unordered_set<size_t> all_up_states(size_t state);
    std::unordered_set<size_t> all_up_states_r(size_t state);
    const std::vector<Edge> & up_edges(size_t state)const {return up_graph[state];};
    const std::vector<Edge> & up_inv_edges(size_t state)const {return up_inv_graph[state];};
    const std::vector<Edge> & down_edges(size_t state)const {return down_graph[state];};
    const std::vector<Edge> & down_inv_edges(size_t state)const {return down_inv_graph[state];};

    CHGraphRT(std::string fname, size_t num_obj=2);


    size_t get_num_of_objectives()const {return num_obj;}
};

class CHShortestPathHeuristic {
private:
    size_t                  source;
    size_t num_obj;

    void compute(size_t cost_idx, const CHGraphRT& chg, const std::unordered_set<size_t>& up_states, const std::unordered_set<size_t>& down_states);
public:

    std::vector<std::unordered_map<size_t, size_t>> data;
  
    CHShortestPathHeuristic(size_t source, const CHGraphRT& chg, const std::unordered_set<size_t>& up_states, const std::unordered_set<size_t>& down_states);
    std::vector<size_t> operator()(size_t node_id);
};

class BOAStarCH {
protected:
    CHGraphRT &chg;


public:
    std::clock_t start_time;
    std::clock_t heuristic_time;
    std::clock_t total_time;
    std::clock_t preprocess_time;
    size_t num_generation = 0;
    size_t num_expansion= 0;
    virtual std::string get_solver_name() {return "BOA* + CH"; }

    BOAStarCH(CHGraphRT &chg):chg(chg){};

    void operator()(size_t source, size_t target, SolutionSet &solutions, unsigned int time_limit=UINT_MAX);

};

class ApexSearchCH {
protected:
    CHGraphRT &chg;
    EPS eps;
    // std::unique_ptr<DominanceChecker> local_dom_checker;
    std::unique_ptr<DominanceChecker> solution_dom_checker;
    MergeStrategy ms=MergeStrategy::SMALLER_G2;

    virtual void insert(ApexPathPairPtr &pp, APQueueHash &queue);
    bool is_dominated(ApexPathPairPtr ap);
    void merge_to_solutions(const ApexPathPairPtr &pp, ApexPathSolutionSet &solutions);

    std::unordered_map<size_t, size_t> min_g2;

public:
    std::clock_t start_time;
    std::clock_t heuristic_time;
    std::clock_t total_time;
    std::clock_t preprocess_time;
    size_t num_generation = 0;
    size_t num_expansion= 0;
    virtual std::string get_solver_name() {return "A*pex + CH"; }

    ApexSearchCH(CHGraphRT &chg, EPS eps):chg(chg), eps(eps){};

    void operator()(size_t source, size_t target, SolutionSet &solutions, unsigned int time_limit=UINT_MAX);

};


