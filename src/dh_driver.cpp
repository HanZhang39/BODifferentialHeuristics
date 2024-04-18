#include <set>
#include <iostream>
#include <algorithm>
#include <vector>
#include <time.h>
#include "boost/range/algorithm/upper_bound.hpp"
#include "MultiVectorHeuristic.h"
#include "ShortestPathHeuristic.h"
#include "Utils/Definitions.h"
#include "Utils/IOUtils.h"
#include "BOAStar.h"

#include <boost/program_options.hpp>
#include<boost/tokenizer.hpp>

int main(int argc, char** argv){
    namespace po = boost::program_options;

    std::vector<std::string> objective_files;

    // Declare the supported options.
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("query,q", po::value<std::string>()->required(), "query file")
        ("map,m",po::value< std::vector<std::string> >(&objective_files)->multitoken(), "files for edge weight")
        ("cutoffTime,t", po::value<int>()->default_value(300), "cutoff time (seconds)")
        ("output,o", po::value<std::string>()->required(), "Name of the output file")
        ("dh", po::value<std::string>()->required(), "Name of the dh")
        ("dhp", po::value<std::string>()->default_value(""), "Name of the dh")
        ("nl", po::value<int>()->default_value(10000000), "Num of landmarks")
        ("update_threshold", po::value<double>()->default_value(0.00001), "Num of landmarks")
        ("update_interval", po::value<size_t>()->default_value(1000000000), "Num of landmarks")
        ("useall", po::value<bool>()->default_value(false), "Using all landmarks")
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);

    if (vm.count("help")) {
        std::cout << desc << std::endl;
        return 1;
    }

    po::notify(vm);

    size_t graph_size;
    std::vector<Edge> edges;

    for (auto file:objective_files){
        std::cout << file << std::endl;
    }


    if (load_gr_files(objective_files, edges, graph_size) == false) {
        std::cout << "Failed to load gr files" << std::endl;
        return -1;
    }

    std::vector<std::pair<size_t, size_t>> queries;
    load_queries(vm["query"].as<std::string>(), queries);

    std::ofstream stats;
    stats.open(vm["output"].as<std::string>(), std::fstream::app);


    AdjacencyMatrix graph(graph_size, edges);
    AdjacencyMatrix inv_graph(graph_size, edges, true);


    std::vector<DHDataPtr> all_data;
    std::vector<DHDataPtr> all_data_path;
    long long cnt = 0;
    for (auto fname :read_filename(vm["dh"].as<std::string>())){
        // std::cout << fname << std::endl;
        all_data.push_back(load_dh(fname));
        for (auto &l: *(all_data.back())  ){
            cnt += l.first.size();
        }
        if (all_data.size() >= vm["nl"].as<int>()){
            break;
        }
    }
    std::cout << "loaded " << all_data.size() << "apex files" << std::endl;
    if ( vm["dhp"].as<std::string>() == ""){
        all_data_path = all_data;
    } else {
        for (auto fname :read_filename(vm["dhp"].as<std::string>())){
            // std::cout << fname << std::endl;
            all_data_path.push_back(load_dh(fname));
            for (auto &l: *(all_data_path.back()) ){
                cnt += l.first.size();
            }

            if (all_data_path.size() >= vm["nl"].as<int>()){
                break;
            }

        }
        std::cout << "loaded " << all_data_path.size() << "path files" << std::endl;
    }
    std::cout << "total size: "<< cnt << std::endl;

    clock_t start_time;
    clock_t runtime;
    SolutionSet sols;

    for (auto entry:queries){
        size_t start = entry.first;
        size_t goal = entry.second;

        stats << start << "\t" << goal << "\t";

        ShortestPathHeuristic sp_heuristic(goal, graph_size, inv_graph);
        Heuristic heuristic = std::bind( &ShortestPathHeuristic::operator(), sp_heuristic, std::placeholders::_1);

        DHBiobjDA H(goal, & sp_heuristic, {all_data}, {all_data_path});
        // DHBiobjDANaiveComax H(goal, & sp_heuristic, {all_data}, {all_data_path});
        H.set_update_interval(vm["update_interval"].as<size_t>() );
        H.set_update_threshold(vm["update_threshold"].as<double>() );

        if ( vm["useall"].as<bool>() ){
            H.init_activation_all();
        } else {
            H.init_activation(start);
        }

        BOAStarMVH solver(graph, {0,0});
        solver.mvh = &H;
        sols.clear();

        start_time =std::clock();
        solver((size_t)start, (size_t)goal, heuristic, sols);
        runtime = std::clock() - start_time;
        std::cout << "Node expansion: " << solver.get_num_expansion() << std::endl;
        std::cout << "Runtime: " <<  ((double) runtime) / CLOCKS_PER_SEC<< std::endl;
        std::cout << "#sols: " <<  sols.size() << std::endl;

        stats << solver.get_num_expansion()
              << "\t" << ((double) runtime) / CLOCKS_PER_SEC
              << "\t" << ((double) H.runtime) / CLOCKS_PER_SEC
              << "\t" << sols.size() << std::endl;

    }

    return 0;
}
