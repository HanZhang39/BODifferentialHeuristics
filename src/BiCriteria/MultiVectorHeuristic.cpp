#include "MultiVectorHeuristic.h"

#include <fstream>
#include <boost/tokenizer.hpp>


std::pair<size_t, size_t> MultiVectorHeuristicBiobj::heuristic(size_t state, std::vector<size_t> g_val, int f2min){
    auto & vec = *get_h_vectors(state);
    if (vec.size() ==0 ) {
        return {0, 0};
    }
    int h2_bound = f2min - g_val[1];
    if ((int)vec[0].second < h2_bound){
        return vec[0];
    }
    if (vec.back().second >= h2_bound){
        return vec.back();
    }

    size_t lb=0, ub = vec.size() - 1;
    while (ub -lb > 1){
        size_t mid = (ub + lb) / 2;
        if (vec[mid].second >= h2_bound){
            lb = mid;
        } else {
            ub = mid;
        }
    }
    return vec[ub];
}


std::vector<std::pair<size_t, size_t>> * DHBiobj::get_h_vectors(size_t state){
    if ( cache.find(state) == cache.end() ){
        auto start_time =std::clock();
        compute_h(state);
        runtime += std::clock() - start_time;
    }
    return cache[state];
}

std::vector<std::pair<size_t, size_t>> * DHBiobjDA::get_h_vectors(size_t state){
    if ( cache.find(state) == cache.end() ){
        auto start_time =std::clock();

        cnt_states += 1;

        if (cnt_states % update_interval == 0){
            init_activation(state, 1);
        }

        compute_h(state);
        runtime += std::clock() - start_time;
    }

    if (cache_timestep[state] != all_data.size()){
        auto start_time =std::clock();
        compute_h(state);
        runtime += std::clock() - start_time;
    }

    return cache[state];
}


std::vector<std::pair<size_t, size_t>> comax(std::vector<std::pair<size_t, size_t>>* vec1,
                                             std::vector<std::pair<size_t, size_t>>* vec2){
    std::vector<std::pair<size_t, size_t>> result;
    size_t prev_y = 0;
    int i = vec1->size() - 1;
    int j = vec2->size() - 1;
    while (i >=0 ||j >= 0){
        if (i<0 || (j>= 0 && vec1->at(i) <= vec2->at(j)) ){
            if ( prev_y < vec2->at(j).second){
                result.push_back(vec2->at(j));
                prev_y = vec2->at(j).second;
            }
            j-= 1;

        } else {
            if ( prev_y < vec1->at(i).second){
                result.push_back(vec1->at(i));
                prev_y = vec1->at(i).second;
            }
            i-= 1;
        }
    }
    std::reverse(result.begin(), result.end());
    return result;
}


std::vector<std::pair<size_t, size_t>> comax_naive(std::vector<std::pair<size_t, size_t>>* vec1,
                                             std::vector<std::pair<size_t, size_t>>* vec2){
    std::vector<std::pair<size_t, size_t>> all_points;


    for (auto & p1: *vec1){
        for (auto & p2: *vec2){
            all_points.push_back({std::max(p1.first, p2.first), std::max(p1.second, p2.second)});
        }
    }

    std::sort(all_points.begin(), all_points.end());
    std::vector<std::pair<size_t, size_t>> result;
    size_t min_g2 = all_points.front().second + 1;

    for (auto & point: all_points){
        if (point.second < min_g2){
            result.push_back(point);
            min_g2 = point.second;
        }
    }

    return result;
}


std::vector<std::pair<size_t, size_t>>* comax(std::vector<std::pair<size_t, size_t>>& vec){
    std::sort(vec.begin(), vec.end());
    auto res = new std::vector<std::pair <size_t, size_t>>();
    for (int i = vec.size() - 1; i >= 0; i--){
        if (res->empty() ||
            (res->back().first > vec[i].first && res->back().second < vec[i].second )
            ){
            res->push_back(vec[i]);
        }
    }
    return res;
}

void DHBiobjDANaiveComax::compute_h(size_t state){
    auto h_val = (*sph)(state);
    size_t min_h1 = h_val[0];
    size_t min_h2 = h_val[1];

    size_t state_2 = goal_state;
    size_t state_1 = state;

    // std::vector<std::vector<std::pair<size_t, size_t>>> all_points;
    std::vector<std::pair<size_t, size_t>> points = {{min_h1, min_h2}};

    size_t i_begin = 0;
    if (cache.find(state) != cache.end()){
        points = * cache[state];
        i_begin = cache_timestep[state];
    }


    if (use_fwd){
    for (size_t i_lm = i_begin; i_lm < all_data.size(); i_lm++){

        auto &data_ptr = all_data[i_lm]->at(state_2);
        auto & data_path_ptr = all_data_path[i_lm]->at(state_1);

        auto max_h1 = data_ptr.first.back() - data_path_ptr.first.front();
        auto max_h2 = data_ptr.second.front() - data_path_ptr.second.back();
        if (max_h1 <= min_h1 || max_h2 <= min_h2){
            continue;
        }
        size_t first_to_check = 0;

        for (size_t j =0; j < data_path_ptr.first.size(); j ++){
            std::vector<std::pair<size_t, size_t>> vector;
            for (size_t i = first_to_check; i < data_ptr.first.size(); i ++ ){
                size_t v1 = std::max( (T)min_h1, data_ptr.first[i] - data_path_ptr.first[j]);
                size_t v2 = std::max( (T)min_h2, data_ptr.second[i] - data_path_ptr.second[j]);

                vector.push_back({v1, v2});
            }
            // do comax
            points = comax_naive(&points, &vector);
        }
    }}

    if (use_bwd){
    std::swap(state_1, state_2);
    for (size_t i_lm = i_begin; i_lm < all_data.size(); i_lm++){
        auto &data_ptr = all_data[i_lm]->at(state_2);
        auto & data_path_ptr = all_data_path[i_lm]->at(state_1);

        auto max_h1 = data_ptr.first.back() - data_path_ptr.first.front();
        auto max_h2 = data_ptr.second.front() - data_path_ptr.second.back();
        if (max_h1 <= min_h1 || max_h2 <= min_h2){
            continue;
        }
        size_t first_to_check = 0;

        for (size_t j =0; j < data_path_ptr.first.size(); j ++){
            std::vector<std::pair<size_t, size_t>> vector;

            for (size_t i = first_to_check; i < data_ptr.first.size(); i ++ ){
                size_t v1 = std::max( (T)min_h1, data_ptr.first[i] - data_path_ptr.first[j]);
                size_t v2 = std::max( (T)min_h2, data_ptr.second[i] - data_path_ptr.second[j]);

                vector.push_back({v1, v2});
            }
            points = comax_naive(&points, &vector);
        }
    }
    }

    if (cache[state] != nullptr){
        delete cache[state];
    }
    
    cache[state] = new std::vector<std::pair<size_t,size_t>>(points.begin(), points.end());
    cache_timestep[state] = all_data.size();


}

void DHBiobj::compute_h(size_t state){
    auto h_val = (*sph)(state);
    size_t min_h1 = h_val[0];
    size_t min_h2 = h_val[1];

    size_t state_2 = goal_state;
    size_t state_1 = state;

    // std::vector<std::vector<std::pair<size_t, size_t>>> all_points;
    std::vector<std::pair<size_t, size_t>> all_points;

    size_t i_begin = 0;
    if (cache.find(state) != cache.end()){
        i_begin = cache_timestep[state];
    }


    if (use_fwd){
    for (size_t i_lm = i_begin; i_lm < all_data.size(); i_lm++){

        auto &data_ptr = all_data[i_lm]->at(state_2);
        auto & data_path_ptr = all_data_path[i_lm]->at(state_1);

        auto max_h1 = data_ptr.first.back() - data_path_ptr.first.front();
        auto max_h2 = data_ptr.second.front() - data_path_ptr.second.back();
        if (max_h1 <= min_h1 || max_h2 <= min_h2){
            continue;
        }
        size_t first_to_check = 0;

        for (size_t j =0; j < data_path_ptr.first.size(); j ++){
            // std::vector<std::pair<size_t, size_t>> diff;
            for (size_t i = first_to_check; i + 1 < data_ptr.first.size(); i ++ ){
                size_t v1 = std::max( (T)min_h1, data_ptr.first[i + 1] - data_path_ptr.first[j]);
                size_t v2 = std::max( (T)min_h2, data_ptr.second[i] - data_path_ptr.second[j]);

                if (v2 == min_h2){break;}
                if (v1 == min_h1){
                    first_to_check = i;
                } else {
                    all_points.push_back({v1, v2});
                }
            }
        }
    }}

    if (use_bwd){
    std::swap(state_1, state_2);
    for (size_t i_lm = i_begin; i_lm < all_data.size(); i_lm++){
        auto &data_ptr = all_data[i_lm]->at(state_2);
        auto & data_path_ptr = all_data_path[i_lm]->at(state_1);

        auto max_h1 = data_ptr.first.back() - data_path_ptr.first.front();
        auto max_h2 = data_ptr.second.front() - data_path_ptr.second.back();
        if (max_h1 <= min_h1 || max_h2 <= min_h2){
            continue;
        }
        size_t first_to_check = 0;

        for (size_t j =0; j < data_path_ptr.first.size(); j ++){
            // std::vector<std::pair<size_t, size_t>> diff;
            for (size_t i = first_to_check; i + 1 < data_ptr.first.size(); i ++ ){
                size_t v1 = std::max((T) min_h1, data_ptr.first[i + 1] - data_path_ptr.first[j]);
                size_t v2 = std::max( (T)min_h2, data_ptr.second[i] - data_path_ptr.second[j]);

                if (v2 == min_h2){break;}
                if (v1 == min_h1){
                    first_to_check = i;
                } else {
                    all_points.push_back({v1, v2});
                }
            }
        }
    }
    }

    if (all_points.size() != 0){
        if (cache.find(state) != cache.end()){
            for (size_t i = 0; i + 1 < cache[state]->size(); i++){
                all_points.push_back({cache[state]->at(i + 1).first, cache[state]->at(i).second});
            }
            delete cache[state];
        }

        auto sorted_vec = comax(all_points);

        std::reverse(sorted_vec->begin(), sorted_vec->end());
        sorted_vec->push_back({0, min_h2});
        for (size_t i =sorted_vec->size() - 1; i >= 1;i --){
            sorted_vec->at(i).first = sorted_vec->at(i-1).first;
        }
        sorted_vec->at(0).first = min_h1;

        cache[state] = sorted_vec;
        cache_timestep[state] = all_data.size();
    } else {
        if (cache.find(state) == cache.end()){
            cache[state] = new std::vector<std::pair<size_t, size_t>>({{min_h1, min_h2}});
        }

        cache_timestep[state] = all_data.size();
    }


}

DHDataPtr load_dh(std::string fname, size_t reserve){
    using namespace boost;
    using namespace std;
    ifstream myfile(fname.c_str());
    if (!myfile.is_open()){
        return nullptr;
    }

    DHDataPtr dh_ptr = make_shared<DHData>();
    dh_ptr->reserve(reserve);

    string line;

    while (getline(myfile, line)){

        if (line[0] == '%'){
            continue;
        }
        // line.erase(std::remove(line.begin(), line.end(), '('), line.end());
        // line.erase(std::remove(line.begin(), line.end(), ')'), line.end());

        tokenizer<char_separator<char>> tok(line, char_separator<char>(";"));
        auto it = tok.begin();

        vector<vector<T>> buffer;
        while (it != tok.end()){
            buffer.emplace_back();
            // cout << "\t" << *it << endl;
            tokenizer<char_separator<char>> tok2(*it, char_separator<char>(","));
            for (auto it2: tok2){
                // cout << "\t\t\t" << atol(it2.c_str()) << endl;
                buffer.back().push_back(atol(it2.c_str()));
            }
            it ++;
        }
        if (buffer.size() == 2){
            dh_ptr->push_back({buffer[0], buffer[1]});
        } else if ( buffer.size() == 0 ){
            dh_ptr->push_back({{}, {}});
        } else {
            cout << "error buffer size " << buffer.size() << endl;
            exit(-1);
        }
    }
    return dh_ptr;
}


std::vector<std::string> read_filename(std::string list_fname){
    using namespace boost;
    using namespace std;
    ifstream myfile(list_fname.c_str());
    if (!myfile.is_open()){
        return {};
    }

    string line;
    tokenizer<char_separator<char>>::iterator beg;
    vector<string> file_list;

    while (getline(myfile, line)){
        file_list.push_back(line);
    }
    return file_list;
}

std::vector<std::pair<size_t, size_t>>* DHBiobjDA::compute_vec(size_t state, size_t i_lm){

    std::vector<std::pair<size_t, size_t>> all_points;
    auto h_val = (*sph)(state);
    size_t min_h1 = h_val[0];
    size_t min_h2 = h_val[1];

    if (use_fwd){
        size_t state_2 = goal_state;
        size_t state_1 = state;
        auto &data_ptr = all_data_pool[i_lm]->at(state_2);
        auto & data_path_ptr = all_data_path_pool[i_lm]->at(state_1);

        auto max_h1 = data_ptr.first.back() - data_path_ptr.first.front();
        auto max_h2 = data_ptr.second.front() - data_path_ptr.second.back();
        if (max_h1 > min_h1 && max_h2 > min_h2){
            size_t first_to_check = 0;
            for (size_t j =0; j < data_path_ptr.first.size(); j ++){
                // std::vector<std::pair<size_t, size_t>> diff;
                for (size_t i = 0; i + 1 < data_ptr.first.size(); i ++ ){
                    size_t v1 = std::max( (T)min_h1, data_ptr.first[i + 1] - data_path_ptr.first[j]);
                    size_t v2 = std::max( (T)min_h2, data_ptr.second[i] - data_path_ptr.second[j]);

                    if (v2 == min_h2){break;}
                    if (v1 == min_h1){
                        first_to_check = i;
                    } else {
                        all_points.push_back({v1, v2});
                    }
                }
            }
        }
    }

    if (use_bwd){
        size_t state_2 = state;
        size_t state_1 = goal_state;

        auto &data_ptr = all_data_pool[i_lm]->at(state_2);
        auto & data_path_ptr = all_data_path_pool[i_lm]->at(state_1);

        auto max_h1 = data_ptr.first.back() - data_path_ptr.first.front();
        auto max_h2 = data_ptr.second.front() - data_path_ptr.second.back();
        if (max_h1 > min_h1 && max_h2 > min_h2){
            size_t first_to_check = 0;

            for (size_t j =0; j < data_path_ptr.first.size(); j ++){
                // std::vector<std::pair<size_t, size_t>> diff;
                for (size_t i = 0; i + 1 < data_ptr.first.size(); i ++ ){
                    size_t v1 = std::max((T) min_h1, data_ptr.first[i + 1] - data_path_ptr.first[j]);
                    size_t v2 = std::max( (T) min_h2, data_ptr.second[i] - data_path_ptr.second[j]);

                    if (v2 == min_h2){break;}
                    if (v1 == min_h1){
                        first_to_check = i;
                    } else {
                        all_points.push_back({v1, v2});
                    }
                }
            }
        }
    }
    auto sorted_vec = comax(all_points);
    std::reverse(sorted_vec->begin(), sorted_vec->end());
    return sorted_vec;
}


double gain(std::vector<size_t>& h_val, std::vector<std::pair<size_t, size_t>>& vec){
    size_t prev = h_val[0];
    size_t offset = h_val[1];
    double sum = 0;
    for (auto& p: vec){
        sum+= (p.first - prev) * (p.second - offset);
        prev = p.first;
    }
    return sum/(h_val[0] * h_val[1]) + 1;
}


void DHBiobjDA::init_activation_all(){
    for (size_t lm = 0; lm < all_data_path_pool.size(); lm++){

        activated[lm] = true;
        all_data.push_back(all_data_pool[lm]);
        all_data_path.push_back(all_data_path_pool[lm]);
    }
}

void DHBiobjDA::init_activation(size_t start, size_t max_to_add){
    std::vector<std::vector<std::pair<size_t, size_t>>*> tmp_vecs(all_data_pool.size(), nullptr);
    for (size_t i = 0; i < tmp_vecs.size(); i ++){
        if (!activated[i]){
            tmp_vecs[i] = compute_vec(start, i);
        }
    }

    auto h_val = (*sph)(start);
    std::vector<std::pair<size_t, size_t>> current_vec;
    if (!all_data.empty()){
        compute_h(start);
        current_vec = *cache[start];
    }
    double curr_val = 1;
    size_t added = 0;
    while (true){
        double max_val = 0;
        int lm = -1;
        for (size_t i = 0; i < tmp_vecs.size(); i++){
            if (activated[i]){
                continue;
            }
            auto merged = comax(&current_vec, tmp_vecs[i]);
            double val_tmp = gain(h_val, merged);
            if (val_tmp > max_val){
                max_val = val_tmp;
                lm = i;
            }
        }


        if (max_val - curr_val >= update_threshold){
            std::cout << "add landmark " << lm << " w/ val: " << max_val << " - " << curr_val << " = " << max_val - curr_val<< std::endl;
             
            activated[lm] = true;
            all_data.push_back(all_data_pool[lm]);
            all_data_path.push_back(all_data_path_pool[lm]);
            current_vec = comax(&current_vec, tmp_vecs[lm]);
            curr_val = max_val;
        } else {
            break;
        }

        if (added >= max_to_add){
            break;
        }
    }

    if (added > 0){
        compute_h(start);
    }

    for (size_t i = 0; i < tmp_vecs.size(); i ++){
        delete tmp_vecs[i];
    }
}
