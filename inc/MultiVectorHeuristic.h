#pragma once

#include "ShortestPathHeuristic.h"
#include <array>

class MultiVectorHeuristicBiobj{
public:
    MultiVectorHeuristicBiobj(){};
    virtual std::pair<size_t, size_t> heuristic(size_t state, std::vector<size_t> g_val, int f2min);
    // return vector of h vals sorted reverse lexicographically
    virtual std::vector<std::pair<size_t, size_t>> * get_h_vectors(size_t state) = 0;



};

typedef int T;
// typedef long long T;
typedef std::vector<std::pair<std::vector<T>,std::vector<T>>> DHData;
typedef std::shared_ptr<DHData> DHDataPtr;


std::vector<std::string> read_filename(std::string list_fname);
DHDataPtr load_dh(std::string fname, size_t reserve=10000);

class DHBiobj: public MultiVectorHeuristicBiobj{
public:
    std::clock_t runtime = 0;

    bool use_fwd = true;
    bool use_bwd = true;

    DHBiobj(size_t goal, ShortestPathHeuristic* sph, std::vector<DHDataPtr> data, std::vector<DHDataPtr> data_path): MultiVectorHeuristicBiobj(),  sph(sph), goal_state(goal), all_data(data), all_data_path(data_path) {};

    virtual ~DHBiobj(){
        for (auto it: cache){
            delete it.second;
        }
    }

    // return the vector of h vals sorted reverse lexicographically
    virtual std::vector<std::pair<size_t, size_t>> * get_h_vectors(size_t state) override;

protected:
    ShortestPathHeuristic* sph;

    size_t goal_state;

    std::unordered_map<size_t, std::vector<std::pair<size_t, size_t>>* > cache;
    std::unordered_map<size_t, size_t> cache_timestep;
    std::vector<DHDataPtr> all_data;
    std::vector<DHDataPtr> all_data_path;


    virtual void compute_h(size_t state);
};


// DA = dynamic activation
class DHBiobjDA: public DHBiobj{
public:

    DHBiobjDA(size_t goal, ShortestPathHeuristic* sph, std::vector<DHDataPtr> data, std::vector<DHDataPtr> data_path):
        DHBiobj(goal, sph, {}, {}),
        activated(data.size(), false), all_data_pool(data),
        all_data_path_pool(data_path)
    {
    };

    // return the vector of h vals sorted reverse lexicographically
    virtual std::vector<std::pair<size_t, size_t>> * get_h_vectors(size_t state) override;

    virtual std::vector<std::pair<size_t, size_t>>* compute_vec(size_t state, size_t lm);

    void init_activation(size_t start, size_t max_to_add=4);
    void init_activation_all();

    void set_update_interval(size_t num){
        update_interval = num;
    }

    void set_update_threshold(double threshold){
        update_threshold = threshold;
    }

protected:

    size_t cnt_states = 0;
    size_t update_interval = 100000000;
    double update_threshold = 0.05;

    std::vector<bool> activated;

    std::vector<DHDataPtr> all_data_pool;
    std::vector<DHDataPtr> all_data_path_pool;

};

// baseline approach of comax
class DHBiobjDANaiveComax: public DHBiobjDA{
public:

    DHBiobjDANaiveComax(size_t goal, ShortestPathHeuristic* sph, std::vector<DHDataPtr> data, std::vector<DHDataPtr> data_path):
        DHBiobjDA(goal, sph, data, data_path)
    {
    };

    // return the vector of h vals sorted reverse lexicographically
    // virtual std::vector<std::pair<size_t, size_t>> * get_h_vectors(size_t state) override;
    
    // 
    // std::vector<std::pair<size_t, size_t>>* compute_vec(size_t state, size_t lm) override;
    
    protected:
    void compute_h(size_t state) override;
};
