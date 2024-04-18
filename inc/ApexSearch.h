#pragma once
#include "Utils/Definitions.h"
#include "Utils/Logger.h"
#include "Utils/MapQueue.h"
#include"DominanceChecker.h"
#include "AbstractSolver.h"


class ApexSearch: public AbstractSolver {
protected:
    size_t num_of_objectives;
    MergeStrategy ms=MergeStrategy::SMALLER_G2;

    std::unique_ptr<DominanceChecker> local_dom_checker;
    std::unique_ptr<DominanceChecker> solution_dom_checker;

    virtual void insert(ApexPathPairPtr &pp, APQueue &queue);
    virtual bool is_dominated(ApexPathPairPtr ap);
    virtual void merge_to_solutions(const ApexPathPairPtr &pp, ApexPathSolutionSet &solutions);
    std::vector<std::vector<ApexPathPairPtr>> expanded;
    void init_search();

    virtual std::string base_alg_name(){return "A*pex";}

public:

    virtual std::string get_solver_name();


    void set_merge_strategy(MergeStrategy new_ms){ms = new_ms;}
    ApexSearch(const AdjacencyMatrix &adj_matrix, EPS eps, const LoggerPtr logger=nullptr);
    virtual void operator()(size_t source, size_t target, Heuristic &heuristic, SolutionSet &solutions, unsigned int time_limit=UINT_MAX) override;
    virtual void operator()(size_t source, size_t target, Heuristic &heuristic, std::vector<ApexPathPairPtr> &solutions, unsigned int time_limit=UINT_MAX);
};



class Dpex: public ApexSearch {
protected:
    void merge_to_solutions(const ApexPathPairPtr &pp, ApexPathSolutionSet &solutions) override ;
    std::string base_alg_name() override{
        return "Dpex";
    }

    bool is_dominated(ApexPathPairPtr ap);

public:


    void set_merge_strategy(MergeStrategy new_ms){ms = new_ms;}
    Dpex(const AdjacencyMatrix &adj_matrix, EPS eps, const LoggerPtr logger=nullptr):ApexSearch(adj_matrix, eps, logger) {};
    virtual void operator()(size_t source, size_t target, Heuristic &heuristic, SolutionSet &solutions, unsigned int time_limit=UINT_MAX) override {
        std::cout << "this function is not implemented for DPex" << std::endl;
    }

    void operator()(size_t source, std::vector<std::vector<std::pair<size_t, size_t>>> &solutions, unsigned int time_limit=UINT_MAX) override;
    void operator()(size_t source, std::unordered_map<size_t, std::vector<ApexPathPairPtr>> &solutions, unsigned int time_limit=UINT_MAX);

};
