#include <Rcpp.h>
#include <map>
#include <string>
#include <vector>
#include <algorithm>
#include <random>

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]  
Rcpp::List get_contigs(const std::vector<std::string> read_kmers, const int dbg_kmer, const int seed){
    // init prefix and suffix vectors
    std::vector<std::string> prefix(read_kmers.size());
    std::vector<std::string> suffix(read_kmers.size());

    for(int i = 0; i < read_kmers.size(); i++){
        // extract single kmer from string vector
        std::string read_kmer = read_kmers[i];

        // extract prefix and suffix of kmer
        prefix[i] = read_kmer.substr(0, dbg_kmer-1);
        suffix[i] = read_kmer.substr(1, dbg_kmer);
    }

    // init hash map of prefix and suffix
    std::map<std::string, std::vector<std::string>> dict;

    // fill hash map with pre- and suffixes
    for(int i = 0; i < prefix.size(); i++){
        std::string key = prefix[i];
        std::string value = suffix[i];

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
    Rcpp::CharacterVector temp_prefix = Rcpp::wrap(prefix);
    Rcpp::CharacterVector temp_suffix = Rcpp::wrap(suffix);
    std::vector<std::string> unique_vals = Rcpp::as<std::vector<std::string>>(
        Rcpp::unique(temp_prefix)
    );
    std::vector<std::string> unique_suffix = Rcpp::as<std::vector<std::string>>(
        Rcpp::unique(temp_suffix)
    );
    unique_vals.insert(unique_vals.end(), unique_suffix.begin(), unique_suffix.end());

    // init balanced counts to zeros
    for(int i = 0; i < unique_vals.size(); i++){
        std::string key = unique_vals[i];
        balanced_count[key] = std::make_pair(0,0);
    }

    // update balance count per node
    // & calls by reference = any changes made, will change its reference
    for(const auto &pair : dict){
        const std::string &node = pair.first;
        const std::vector<std::string> &edges = pair.second;

        // out-degree
        balanced_count[node].second += edges.size(); 

        // in-degree
        for(const auto &edge : edges){
            const std::string &current_edge = edge;
            balanced_count[current_edge].first += 1;
        }
    }

    // find branching points
    std::vector<std::string> keys_with_balanced_count;
    for(const auto &pair : balanced_count){
        if(pair.second.first != 1 || pair.second.second != 1){
            // keep branching point only if it's a node
            if(dict.count(pair.first)){
                keys_with_balanced_count.push_back(pair.first);
            }
        }
    }

    // get contigs
    std::vector<std::string> contigs;
    for(const auto &node : keys_with_balanced_count){
        const std::vector<std::string> edges = dict[node];
        for(const auto &edge : edges){
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
    std::vector<std::string> unique_contigs = Rcpp::as<std::vector<std::string>>(
        Rcpp::unique(r_contigs)
    );

    // Randomly shuffle the order of contigs for assembly
    std::mt19937 engine(seed);
    Rcpp::List contig_matrix(10);
    for(int i = 0; i < contig_matrix.size(); i++){
        std::vector<std::string> contigs_copy = unique_contigs;
        std::shuffle(contigs_copy.begin(), contigs_copy.end(), engine);
        contig_matrix[i] = contigs_copy;
    }
    return contig_matrix;
}

// [[Rcpp::export]]
std::vector<std::string> assemble_contigs(Rcpp::List contig_matrix, const int dbg_kmer){
    // Loop over all the contig subsets and return assemblies in-place.
    // Note. Below code looks worse than it appears. In reality, it executes fast
    for(int contig_ind = 0; contig_ind < contig_matrix.size(); contig_ind++){
        Rcpp::StringVector contigs = contig_matrix[contig_ind];
        for(int kmer = (dbg_kmer-1); kmer > 0; kmer--){
            bool len_changed = true;
            while(len_changed){
                int temp = contigs.size();
                for(int i = 0; i < contigs.size(); i++){
                    if (contigs[i] == "") continue;
                    for(int j = contigs.size()-1; j >= 0; j--){
                        if(contigs[i] != contigs[j]){
                            std::string str_contig_i = Rcpp::as<std::string>(contigs[i]);
                            std::string str_contig_j = Rcpp::as<std::string>(contigs[j]);
                            std::string suffix = str_contig_i.substr(
                                contigs[i].size()-kmer, 
                                contigs[i].size()
                            );
                            std::string prefix = str_contig_j.substr(0, kmer);
                            if(suffix == prefix){
                                std::string no_prefix = str_contig_j.substr(
                                    kmer, 
                                    contigs[j].size()
                                );
                                contigs[i] = str_contig_i.append(no_prefix);
                                contigs[j] = "";
                            }
                        }
                    }
                }
                for(int i = contigs.size()-1; i >= 0; i--){
                    if(contigs[i] == ""){
                        contigs.erase(i);
                    }
                }
                len_changed = temp != contigs.size();
            }
        }
        contig_matrix[contig_ind] = contigs;
    }

    // flatten list into character vector and discard duplicates
    std::vector<std::string> flat_contig_matrix;
    for(int i = 0; i < contig_matrix.size(); i++){
        std::vector<std::string> current_contigs = contig_matrix[i];
        for(int j = 0; j < current_contigs.size(); j++){
            flat_contig_matrix.push_back(current_contigs[j]);
        }
    }

    // Discard duplicates
    Rcpp::CharacterVector r_flat_contig_matrix = Rcpp::wrap(flat_contig_matrix);
    std::vector<std::string> contig_matrix_set = Rcpp::as<std::vector<std::string>>(
        Rcpp::unique(r_flat_contig_matrix)
    );
    return contig_matrix_set;
}

/*
    Get reverse complement of a DNA string
*/
std::string reverse_complement(const std::string &sequence) {
    // reverse complement map
    std::unordered_map<char, char> complement = {
        {'A', 'T'}, 
        {'C', 'G'}, 
        {'G', 'C'}, 
        {'T', 'A'}
    };
    std::string reverse_seq(sequence.rbegin(), sequence.rend());
    for(char &c : reverse_seq){
        c = complement[c];
    }
    return reverse_seq;
}

// [[Rcpp::export]]
Rcpp::List calc_breakscore(std::vector<std::string> path, std::vector<std::string> sequencing_reads, const int kmer, Rcpp::DataFrame bp_table){
    // extract columns from Rcpp::DataFrame
    std::vector<std::string> bp_kmer = Rcpp::as<std::vector<std::string>>(bp_table["kmer"]);
    std::vector<double> bp_prob = Rcpp::as<std::vector<double>>(bp_table["prob"]);
    
    // hash map of kmers and probability values
    std::map<std::string, double> bp_matrix;
    for(int i = 0; i < bp_kmer.size(); i++){
        bp_matrix[bp_kmer[i]] = bp_prob[i];
    }

    // init break score vector
    std::vector<double> break_score_vector(path.size());

    for(int i = 0; i < path.size(); i++){
        // init break_score
        double break_score = 0;
        std::string current_path = path[i];

        for(const auto &read : sequencing_reads){
            // find exact match
            std::size_t pos = current_path.find(read);

            if(pos != std::string::npos){
                // get start pos of broken kmer
                // take max of start pos of kmer or start of path
                int start_pos_ind = std::max(0, (int)pos-(kmer/2));

                // extract broken kmer
                std::string broken_kmer = current_path.substr(start_pos_ind, kmer);

                // get corresponding probability from hash map
                double prob = bp_matrix[broken_kmer];

                // get reverse complement kmer if no match found
                if(prob == 0){
                    std::string rev_comp_kmer = reverse_complement(broken_kmer);
                    prob = bp_matrix[rev_comp_kmer];
                }

                // update probability score
                break_score += prob;
            }
        }
        break_score_vector[i] = break_score;
    }

    return Rcpp::List::create(
        Rcpp::Named("sequence") = path,
        Rcpp::Named("bp_score") = break_score_vector
    );
}