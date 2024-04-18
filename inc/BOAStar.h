#ifndef BI_CRITERIA_BOA_STAR_H
#define BI_CRITERIA_BOA_STAR_H

#include <vector>
#include "Utils/Definitions.h"
#include "Utils/Logger.h"
#include "AbstractSolver.h"
#include "MultiVectorHeuristic.h"


class BOAStar: public AbstractSolver {
protected:
    std::clock_t start_time;

    std::vector<std::pair<std::clock_t, NodePtr>> solution_log;
    void log_solution(NodePtr);

public:
    virtual std::string get_solver_name() {return "BOA*"; }

    BOAStar(const AdjacencyMatrix &adj_matrix, Pair<double> eps, const LoggerPtr logger=nullptr);

    void operator()(size_t source, size_t target, Heuristic &heuristic, SolutionSet &solutions, unsigned int time_limit=UINT_MAX);

    std::vector<std::pair<std::clock_t, NodePtr>> get_sol_log(){return solution_log;}
};

class NAMOA: public BOAStar {
public:

    NAMOA(const AdjacencyMatrix &adj_matrix, Pair<double> eps, const LoggerPtr logger=nullptr):BOAStar(adj_matrix, eps, logger) {};

    void operator()(size_t source, size_t target, Heuristic &heuristic, SolutionSet &solutions, unsigned int time_limit=UINT_MAX) override;

};

class BOAStarMVH: public BOAStar {
public:
    // std::unique_ptr<MultiVectorHeuristicBiobj> mvh;
    MultiVectorHeuristicBiobj* mvh = nullptr;

    BOAStarMVH(const AdjacencyMatrix &adj_matrix, Pair<double> eps, const LoggerPtr logger=nullptr):BOAStar(adj_matrix, eps, logger) {};

    void operator()(size_t source, size_t target, Heuristic &heuristic, SolutionSet &solutions, unsigned int time_limit=UINT_MAX) override;

};

class BOD: public BOAStar {
protected:
public:
    std::string get_solver_name() override{
        return "BOD";
    }
    BOD(const AdjacencyMatrix &adj_matrix,  const LoggerPtr logger=nullptr): BOAStar(adj_matrix, {0,0}, logger) {};

    virtual void operator()(size_t source, size_t target, Heuristic &heuristic, SolutionSet &solutions, unsigned int time_limit=UINT_MAX) override {
        std::cout << "this function is not implemented for DPex" << std::endl;
    }
    void operator()(size_t source, std::vector<std::vector<std::pair<size_t, size_t>>> &solutions, unsigned int time_limit=UINT_MAX) override;

};



#endif //BI_CRITERIA_BOA_STAR_H
