#include <iostream>
#include <memory>
#include <time.h>
#include <fstream>

#include "ShortestPathHeuristic.h"
#include "Utils/Definitions.h"
#include "Utils/IOUtils.h"
#include "Utils/Logger.h"
#include "BOAStar.h"
#include "PPA.h"
#include "SingleCriteria.h"
#include "ApexSearch.h"
#include "NAMOA.h"

#include <boost/program_options.hpp>
#include<boost/tokenizer.hpp>

using namespace std;


const MergeStrategy DEFAULT_MERGE_STRATEGY = MergeStrategy::SMALLER_G2;


void single_run_map(size_t graph_size, std::vector<Edge> & edges, size_t source, std::string output_file, std::string algorithm, LoggerPtr logger, MergeStrategy merge_strategy, double eps, int time_limit) {

    AdjacencyMatrix graph(graph_size, edges);

    std::unique_ptr<AbstractSolver> solver;
    if (algorithm == "BOD"){
        solver = std::make_unique<BOD>(graph, logger);
    }else if (algorithm == "Dpex"){
        EPS eps_vec (graph.get_num_of_objectives(), eps);
        solver = std::make_unique<Dpex>(graph, eps_vec, logger);
        ((Dpex*) solver.get())->set_merge_strategy(merge_strategy);
    }else{
        std::cerr << "unknown solver name" << std::endl;
        exit(-1);
    }

    vector<vector<std::pair<size_t, size_t>>> solutions;

    auto start =std::clock();
    (*solver)(source, solutions, time_limit);
    auto runtime = std::clock() - start;

    int num_exp, num_gen;
    num_exp = solver->get_num_expansion();
    num_gen = solver->get_num_generation();

    std::cout << "Node expansion: " << num_exp << std::endl;
    std::cout << "Runtime: " <<  ((double) runtime) / CLOCKS_PER_SEC<< std::endl;

    // output to file
    std::ofstream stats_fstream;
    stats_fstream.open(output_file, std::fstream::out);

    stats_fstream << "% " << algorithm << "-" << algorithm << " (" << eps << ")" << "\t"
                  << source << "\t"  << "\t"
                  << num_gen << "\t"
                  << num_exp << "\t"
                  << solutions.size() << "\t"
                  << (double) runtime / CLOCKS_PER_SEC << "\t"
                  << std::endl;
    for (auto& l : solutions){
        bool need_to_add_comma = false;
        for (auto node: l){
            if (need_to_add_comma){stats_fstream<< ", ";}
            stats_fstream << node.first ;
            need_to_add_comma = true;
        }
        stats_fstream << ";";
        need_to_add_comma = false;
        for (auto node: l){
            if (need_to_add_comma){stats_fstream<< ", ";}
            stats_fstream << node.second;
            need_to_add_comma = true;
        }
        stats_fstream << endl;
    }
}

int main(int argc, char** argv){
    namespace po = boost::program_options;

    std::vector<string> edge_cost_files;

    // Declare the supported options.
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("start,s", po::value<int>()->default_value(-1), "start location")
        ("map,m",po::value< std::vector<string> >(&edge_cost_files)->multitoken(), "files for edge costs")
        ("eps,e", po::value<double>()->default_value(0), "approximation factor")
        ("algorithm,a", po::value<std::string>()->default_value("Apex"), "solvers ()")
        ("merge", po::value<std::string>()->default_value("MORE_SLACK"), "strategy for merging apex node pair: SMALLER_G2, RANDOM or MORE_SLACK")
        ("cutoffTime,t", po::value<int>()->default_value(300), "cutoff time (seconds)")
        ("output,o", po::value<std::string>()->required(), "Name of the output file")
        // ("logging_file", po::value<std::string>()->default_value(""), "logging file" )
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);

    if (vm.count("help")) {
        std::cout << desc << std::endl;
        return 1;
    }

    po::notify(vm);
    srand((int)time(0));

    LoggerPtr logger = nullptr;

    // Load files
    size_t graph_size;
    std::vector<Edge> edges;

    for (auto file:edge_cost_files){
        std::cout << file << std::endl;
    }

    if (load_gr_files(edge_cost_files, edges, graph_size) == false) {
        std::cout << "Failed to load edge cost files" << std::endl;
        return -1;
    }

    std::cout << "Graph Size: " << graph_size << std::endl;

    MergeStrategy merge_strategy = DEFAULT_MERGE_STRATEGY;

    if(vm["merge"].as<std::string>() == "SMALLER_G2"){
        merge_strategy = MergeStrategy::SMALLER_G2;
    }else if(vm["merge"].as<std::string>() == "SMALLER_G2_FIRST"){
        merge_strategy = MergeStrategy::SMALLER_G2_FIRST;
    }else if(vm["merge"].as<std::string>() == "RANDOM"){
        merge_strategy = MergeStrategy::RANDOM;
    }else if(vm["merge"].as<std::string>() == "MORE_SLACK"){
        merge_strategy = MergeStrategy::MORE_SLACK;
    }else if(vm["merge"].as<std::string>() == "REVERSE_LEX"){
        merge_strategy = MergeStrategy::REVERSE_LEX;
    }else{
        std::cerr << "unknown merge strategy" << std::endl;
    }

    single_run_map(graph_size, edges, vm["start"].as<int>(), vm["output"].as<std::string>(), vm["algorithm"].as<std::string>(), logger, merge_strategy, vm["eps"].as<double>(), vm["cutoffTime"].as<int>());

    return 0;
}
