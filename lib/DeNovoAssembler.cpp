#include <Rcpp.h>
#include <map>
#include <string>
#include <vector>
#include <algorithm>

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]  
Rcpp::List get_prefix_and_suffix(Rcpp::CharacterVector read_kmers, int dbg_kmer){
    // init prefix and suffix vectors
    Rcpp::CharacterVector prefix(read_kmers.size());
    Rcpp::CharacterVector suffix(read_kmers.size());

    for(int i = 0; i < read_kmers.size(); i++){
        // extract single kmer from string vector
        std::string read_kmer = Rcpp::as<std::string>(read_kmers[i]);

        // extract prefix and suffix of kmer
        prefix[i] = read_kmer.substr(0, dbg_kmer-1);
        suffix[i] = read_kmer.substr(1, dbg_kmer);
    }
    return Rcpp::List::create(Rcpp::Named("prefix") = prefix,
                              Rcpp::Named("suffix") = suffix); 
}

// [[Rcpp::export]]  
Rcpp::CharacterVector get_contigs(const Rcpp::List prefix_suffix){
    Rcpp::CharacterVector prefix = prefix_suffix["prefix"];
    Rcpp::CharacterVector suffix = prefix_suffix["suffix"];

    // init hash map of prefix and suffix
    std::map<std::string, std::vector<std::string>> dict;

    // fill hash map with pre- and suffixes
    for(int i = 0; i < prefix.size(); i++){
        std::string key = Rcpp::as<std::string>(prefix[i]);
        std::string value = Rcpp::as<std::string>(suffix[i]);

        if(dict[key].empty()){
            // creates new entry in map if previously not existing
            dict[key].push_back(value);
        } else {
            // points to iterator if found, else end of vector
            auto it = std::find(dict[key].begin(), dict[key].end(), value);
            if(it == dict[key].end()){
                // if value not present
                dict[key].push_back(value);
            }
        }
    }

    // init hash map of balanced counts
    std::map<std::string, std::pair<int, int>> balanced_count;

    // get unique prefix and suffix
    std::vector<std::string> unique_vals = Rcpp::as<std::vector<std::string>>(Rcpp::unique(prefix));
    std::vector<std::string> unique_suffix = Rcpp::as<std::vector<std::string>>(Rcpp::unique(suffix));
    unique_vals.insert(unique_vals.end(), unique_suffix.begin(), unique_suffix.end());

    // init balanced counts to zeros
    for(int i = 0; i < unique_vals.size(); i++){
        std::string key = unique_vals[i];
        balanced_count[key] = std::make_pair(0,0);
    }

    // update balance count per node
    for(auto const& pair : dict){
        std::string node = pair.first;
        auto edges = pair.second;

        // out-degree
        balanced_count[node].second += edges.size(); 

        // in-degree
        for(auto const& edge : edges){
            std::string current_edge = edge;
            balanced_count[current_edge].first += 1;
        }
    }

    // find branching points
    std::vector<std::string> keys_with_balanced_count;
    for(auto const& pair : balanced_count){
        if(pair.second.first != 1 || pair.second.second != 1){

            // keep branching point only if it's a node
            if(dict.count(pair.first)){
                keys_with_balanced_count.push_back(pair.first);
            }
        }
    }

    // get contigs
    std::vector<std::string> contigs;
    for(auto const& node : keys_with_balanced_count){
        std::vector<std::string> edges = dict[node];
        for(auto const& edge : edges){
            std::string current_node = edge;
            std::string path = node;       

            while(std::find(keys_with_balanced_count.begin(), 
                            keys_with_balanced_count.end(), 
                            current_node) == keys_with_balanced_count.end()){
                if(dict[current_node].empty()) break;
                path.append(current_node.substr(current_node.size()-1, 1));
                current_node = dict[current_node][0];
            }
            path.append(current_node.substr(current_node.size()-1, 1));
            contigs.push_back(path);
        }
    }

    // To use Rcpp::unique, need to convert to Rcpp Vector types
    Rcpp::CharacterVector r_contigs = Rcpp::wrap(contigs);
    Rcpp::CharacterVector unique_contigs = Rcpp::unique(r_contigs);
    return unique_contigs;
}

// [[Rcpp::export]]
Rcpp::List assemble_contigs(List contig_matrix, int dbg_kmer){

    // Loop over all the contig subsets
    for(int contig_ind = 0; contig_ind < contig_matrix.size(); contig_ind++){
        Rcpp::StringVector contigs = contig_matrix[contig_ind];
        for(int kmer = (dbg_kmer-1); kmer > 0; kmer--){
            bool len_changed = true;
            while(len_changed){
                int temp = contigs.size();
                for(int i = 0; i < contigs.size(); i++){
                    if (contigs[i] == "") continue;
                    for(int j = contigs.size()-1; j >= 0; j--){
                        if (contigs[i] != contigs[j]) {
                            std::string str_contig_i = Rcpp::as<std::string>(contigs[i]);
                            std::string str_contig_j = Rcpp::as<std::string>(contigs[j]);
                            std::string suffix = str_contig_i.substr(contigs[i].size()-kmer, contigs[i].size());
                            std::string prefix = str_contig_j.substr(0, kmer);
                            if (suffix == prefix) {
                                std::string no_prefix = str_contig_j.substr(kmer, contigs[j].size());
                                contigs[i] = str_contig_i.append(no_prefix);
                                contigs[j] = "";
                            }
                        }
                    }
                }
                for(int i = contigs.size()-1; i >= 0; i--){
                    if (contigs[i] == "") {
                        contigs.erase(i);
                    }
                }
                len_changed = temp != contigs.size();
            }
        }
        contig_matrix[contig_ind] = contigs;
    }
    return(contig_matrix);
}