#include "MultiVectorHeuristic.h"
#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <iostream>
#include <fstream>

using namespace std;

inline bool is_bounded(T apex1, T apex2, T cost1, T cost2, double eps){
    return ((double) cost1 <= (1 + eps) * apex1) && ((double) cost2 <= (1 + eps) * apex2);
}

struct BiObjectiveAP {
    T apex1;
    T apex2;
    T cost1;
    T cost2;
    BiObjectiveAP(T cost1, T cost2): apex1(cost1), apex2(cost2), cost1(cost1), cost2(cost2) {};

    void update_apex(T apex1, T apex2){
        this->apex1 = min(this->apex1, apex1);
        this->apex2 = min(this->apex2, apex2);
    }
    void update_cost(T cost1, T cost2){
        this->cost1 = cost1;
        this->cost2 = cost2;
        this->apex1 = min(this->apex1, cost1);
        this->apex2 = min(this->apex2, cost2);
    }
};

std::string format_apex_vec(const std::vector<BiObjectiveAP> &ap_vec) {
    std::stringstream result_stream;
    bool is_first = true;
    for (auto & ap: ap_vec){
        result_stream << (is_first ? "":",") << ap.apex1;
        is_first = false;
    }
    is_first = true;
    result_stream << ";";
    for (auto & ap: ap_vec){
        result_stream << (is_first ? "":",") << ap.apex2;
        is_first = false;
    }

    return result_stream.str();
}

std::string format_cost_vec(const std::vector<BiObjectiveAP> &ap_vec) {
    std::stringstream result_stream;
    bool is_first = true;
    for (auto & ap: ap_vec){
        result_stream << (is_first ? "":",") << ap.cost1;
        is_first = false;
    }
    is_first = true;
    result_stream << ";";
    for (auto & ap: ap_vec){
        result_stream << (is_first ? "":",") << ap.cost2;
        is_first = false;
    }

    return result_stream.str();
}

inline void merge_to_ap_vec(std::vector<BiObjectiveAP>& ap_vec, T& cost1, T& cost2, double & eps){
    if (ap_vec.size() == 0){
        ap_vec.emplace_back(cost1, cost2);
    } else if (is_bounded(ap_vec.back().apex1, ap_vec.back().apex2, cost1, cost2, eps)){
        ap_vec.back().update_cost(cost1, cost2);
    } else if (is_bounded(cost1, cost2, ap_vec.back().cost1, ap_vec.back().cost2, eps)){
        ap_vec.back().update_apex(cost1, cost2);
    } else {
        ap_vec.emplace_back(cost1, cost2);
    }
}

std::vector<BiObjectiveAP> merge_points(std::pair<std::vector<T>,std::vector<T>> & input, double eps){
    if (input.first.size() != input.second.size()){
        std::cout<< "wrong sizes"<< std::endl;
        exit(1);
    }

    std::vector<BiObjectiveAP> result_vec;

    for (size_t i = 0; i < input.first.size(); i ++){
        merge_to_ap_vec(result_vec, input.first[i], input.second[i], eps);
    }

    return result_vec;
}

int main(int argc, char** argv){
    namespace po = boost::program_options;

    // Declare the supported options.
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("input,i", po::value<std::string>()->required(), "start location")
        ("eps,e", po::value<double>()->default_value(0), "approximation factor")
        ("outputApex,o", po::value<std::string>()->required(), "output file for apex")
        ("outputPath,p", po::value<std::string>()->required(), "output file for path")
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);

    if (vm.count("help")) {
        std::cout << desc << std::endl;
        return 1;
    }

    po::notify(vm);

    auto input_data = load_dh(vm["input"].as<std::string>());

    std::ofstream apex_fstream(vm["outputApex"].as<std::string>());
    std::ofstream pathcost_fstream(vm["outputPath"].as<std::string>());

    for(size_t i = 0; i < input_data->size(); i ++){
        auto res = merge_points(input_data->at(i), vm["eps"].as<double>());

        apex_fstream << format_apex_vec(res) << std::endl;
        pathcost_fstream << format_cost_vec(res) << std::endl;
    }

}
