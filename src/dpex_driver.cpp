#include <iostream>
#include <memory>
#include <time.h>
#include <fstream>
#include <boost/program_options.hpp>
#include<boost/tokenizer.hpp>

#include "Utils/IOUtils.h"
#include "ApexSearch.h"

using namespace std;

const MergeStrategy DEFAULT_MERGE_STRATEGY = MergeStrategy::SMALLER_G2;
std::string alg_variant = "";

void run_query(size_t graph_size, std::vector<Edge> & edges, size_t source, vector<size_t> locations, string output_file, MergeStrategy ms, double eps, int time_limit) {
    auto runtime = std::clock();

    AdjacencyMatrix graph(graph_size, edges);
    AdjacencyMatrix inv_graph(graph_size, edges, true);
    cout << eps << endl;
    EPS eps_vec (graph.get_num_of_objectives(), eps);
    LoggerPtr logger = nullptr;

    Dpex solver(graph, eps_vec, logger);
    solver.set_merge_strategy(ms);

    unordered_map<size_t, vector<ApexPathPairPtr>> sols;
    for (auto loc:locations){
        sols[loc] = {};
    }
    auto start =std::clock();
    solver(source, sols, time_limit);
    runtime = std::clock() - start;

    std::cout << "Node expansion: " << solver.get_num_expansion() << std::endl;
    std::cout << "Runtime: " <<  ((double) runtime) / CLOCKS_PER_SEC<< std::endl;
    int num_exp, num_gen;
    num_exp = solver.get_num_expansion();
    num_gen = solver.get_num_generation();

    std::ofstream stats;
    stats.open("log.txt", std::fstream::app);
    stats << "Dpex+" << " (" << eps << ")" << "\t"
           << source << "\t"  << "\t"
           << num_gen << "\t"
           << num_exp << "\t"
           << (double) runtime / CLOCKS_PER_SEC << "\t"
        // << (double) prep_time / CLOCKS_PER_SEC
           << std::endl;


    std::ofstream output;
    output.open(output_file, std::fstream::app);
    for (auto& l : sols){
        size_t k = l.first;

        output << source << "_" << k << ": (";
        bool prev = false;
        for (auto node: l.second){
            if (prev){stats<< "\t";}
            output << node->apex->g << ", ";
            prev = true;
        }
        output << ")" << std::endl;
    }

    std::cout << "-----End Single Example-----" << std::endl;
}

bool load_locations(std::string query_file, std::vector<size_t> &queries_out) {
    std::ifstream   file(query_file.c_str());

    if (file.is_open() == false) {
        return false;
    }

    std::string line;
    while (file.eof() == false) {
        std::getline(file, line);

        if (line == "") {
            break;
        } else if (line[0] == '#') {
            continue; // Commented out queries
        }


        queries_out.push_back(std::stoul(line));
    }
    return true;
}

int main(int argc, char** argv){
    namespace po = boost::program_options;

    std::vector<string> objective_files;

    // Declare the supported options.
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("start,s", po::value<int>()->default_value(-1), "start location")
        ("goals,g", po::value<std::string>()->default_value(""), "goal locations")
        ("map,m",po::value< std::vector<string> >(&objective_files)->multitoken(), "files for edge weight")
        ("eps,e", po::value<double>()->default_value(0), "approximation factor")
        ("merge", po::value<std::string>()->default_value(""), "strategy for merging apex node pair: SMALLER_G2, RANDOM or MORE_SLACK")
        ("cutoffTime,t", po::value<int>()->default_value(300), "cutoff time (seconds)")
        ("output,o", po::value<std::string>()->required(), "Name of the output file")
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);

    if (vm.count("help")) {
        std::cout << desc << std::endl;
        return 1;
    }

    po::notify(vm);
    srand((int)time(0));

    // Load files
    size_t graph_size;
    std::vector<Edge> edges;

    for (auto file:objective_files){
        std::cout << file << std::endl;
    }


    if (load_gr_files(objective_files, edges, graph_size) == false) {
        std::cout << "Failed to load gr files" << std::endl;
        return -1;
    }

    std::cout << "Graph Size: " << graph_size << std::endl;

    // Build graphs
    MergeStrategy ms = DEFAULT_MERGE_STRATEGY;
    alg_variant = vm["merge"].as<std::string>();

    if(vm["merge"].as<std::string>() == "SMALLER_G2"){
        ms = MergeStrategy::SMALLER_G2;
    }else if(vm["merge"].as<std::string>() == "SMALLER_G2_FIRST"){
        ms = MergeStrategy::SMALLER_G2_FIRST;
    }else if(vm["merge"].as<std::string>() == "RANDOM"){
        ms = MergeStrategy::RANDOM;
    }else if(vm["merge"].as<std::string>() == "MORE_SLACK"){
        ms = MergeStrategy::MORE_SLACK;
    }else if(vm["merge"].as<std::string>() == "REVERSE_LEX"){
        ms = MergeStrategy::REVERSE_LEX;
    }else{
        std::cerr << "unknown merge strategy" << std::endl;
    }

    vector<size_t> locs;
    load_locations(vm["goals"].as<std::string>(), locs);
    std::cout << "#locs: " << locs.size() << std::endl;

    run_query(graph_size, edges, vm["start"].as<int>(), locs, vm["output"].as<std::string>(), ms, vm["eps"].as<double>(), vm["cutoffTime"].as<int>());

    return 0;
}
