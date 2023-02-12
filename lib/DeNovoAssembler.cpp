// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins("cpp11")]]
#include <Rcpp.h>
#include <map>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include <random>
#include <cstdio>

// edlib dependency
#include "../lib/edlib/edlib.h"

/**
 * Get reverse complement of a DNA string.
 * @param sequence kmer broken after aligning read with de novo sequence/
 * @return reverse_seq reverse complement of sequence.
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

/**
 * Fast levenshtein distance calculation with edlib.
 * @param query de novo sequence.
 * @param target true solution.
 * @return levenshtein distance between query and target.
*/
int calc_levenshtein(const std::string query, const std::string target){
  // edlib function
  EdlibAlignResult result = edlibAlign(
    query.c_str(), query.length(), 
    target.c_str(), target.length(), 
    edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0)
  );
  
  // return edit distance
  if(result.status == EDLIB_STATUS_OK){
    return result.editDistance;
  }

  return 0;
}

/**
 * Generate contigs from sonicated sequencing reads. Logic flow is:
    * 1) Generate prefix and suffix.
    * 2) Find balanced counts for in- and out-degrees per node.
    * 3) Find branching points.
    * 4) Generate contigs.
    * 5) Shuffle order of contigs and assemble each run.
 * @param read_kmer kmers generated from sequencing reads.
 * @param dbg_kmer size of kmer per de bruijn graph node.
 * @param seed set seed for reproducible results.
 * @return contig_matrix list of shuffled contigs.
*/
// [[Rcpp::export]] 
Rcpp::List get_contigs(const std::vector<std::string> read_kmers, 
                       const int dbg_kmer, 
                       const int seed){
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
    Rcpp::List contig_matrix(100000);
    for(int i = 0; i < contig_matrix.size(); i++){
        std::vector<std::string> contigs_copy = unique_contigs;
        std::shuffle(contigs_copy.begin(), contigs_copy.end(), engine);
        contig_matrix[i] = contigs_copy;
    }
    return contig_matrix;
}

/**
 * Assemble each batch of contigs through brute force alignment.
 * @param contig_matrix list of shuffled contigs.
 * @param dbg_kmer size of kmer per de bruijn graph node.
 * @return top_5_percent_matrix top 5% of assembled solutions by length.
*/
// [[Rcpp::export]]
std::vector<std::string> assemble_contigs(Rcpp::List contig_matrix, 
                                          const int dbg_kmer){
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

    // sort vector in descending order by length
    std::sort(
        contig_matrix_set.begin(),
        contig_matrix_set.end(),
        [](const std::string &seq_a, const std::string &seq_b) -> bool {
            return seq_a.length() > seq_b.length();
    });

    // Retain top 5% of assembled solution by length
    int top_5_percent = contig_matrix_set.size()*0.05;
    std::vector<std::string> top_5_percent_matrix(
        contig_matrix_set.begin(),
        contig_matrix_set.begin()+top_5_percent
    );

    return top_5_percent_matrix;
}

/**
 * Calculate breakage score for each de novo assembled solutions.
 * @param path top 5% of assembled solutions by length.
 * @param sequencing_reads all sequencing reads from ultrasonication.
 * @param true_solution the true solution for this experiment. 
 * @param kmer size of k-meric break.
 * @param bp_table table of kmer and associated breakage probabilities.
 * @return List of de novo assembled solutions and various calculations performed.
*/
// [[Rcpp::export]]
Rcpp::List calc_breakscore(std::vector<std::string> path, 
                           std::vector<std::string> sequencing_reads, 
                           std::string true_solution,
                           const int kmer, 
                           Rcpp::DataFrame bp_table){
    // extract columns from Rcpp::DataFrame
    std::vector<std::string> bp_kmer = Rcpp::as<std::vector<std::string>>(bp_table["kmer"]);
    std::vector<double> bp_prob = Rcpp::as<std::vector<double>>(bp_table["prob"]);

    // hash map of kmers, probability values and break counts
    std::map<std::string, std::pair<double, int>> bp_matrix;
    for(int i = 0; i < bp_kmer.size(); i++){
        bp_matrix[bp_kmer[i]] = std::make_pair(bp_prob[i], 0);
    }

    // init break score vector
    std::vector<double> break_score_vector(path.size());
    std::vector<double> norm_break_score_vector(path.size());
    std::vector<double> break_score_norm_by_len_vector(path.size());
    std::vector<double> nr_of_breaks_vector(path.size());

    for(int i = 0; i < path.size(); i++){
        // init total count of kmer breaks
        int total_breaks = 0;

        for(const auto &read : sequencing_reads){
            // find exact match
            std::size_t pos = path[i].find(read);

            if(pos != std::string::npos){
                // get start pos of broken kmer
                // take max of start pos of kmer or start of path
                int start_pos_ind = std::max(0, (int)pos-(kmer/2));

                // extract broken kmer
                std::string broken_kmer = path[i].substr(start_pos_ind, kmer);

                // Check if the broken kmer is in the bp_matrix
                auto it = bp_matrix.find(broken_kmer);
                if(it != bp_matrix.end()){
                    std::get<1>(it->second)++;
                    total_breaks += 1;
                } else {
                    // Get the reverse complement of the broken kmer
                    std::string rev_comp_kmer = reverse_complement(broken_kmer);

                    // Check if the reverse complement kmer is in the bp_matrix
                    auto it_rev = bp_matrix.find(rev_comp_kmer);
                    if(it_rev != bp_matrix.end()){
                        std::get<1>(it_rev->second)++;
                        total_breaks += 1;
                    }
                }
            }
        }

        // Loop over the bp_matrix
        for(auto it = bp_matrix.begin(); it != bp_matrix.end(); it++){
            double prob = std::get<0>(it->second);
            int break_count = std::get<1>(it->second);

            // Multiply non-zero break_scores with the probability values,
            // add the result to the break_score_vector
            if(break_count != 0){
                double break_score = prob*break_count;
                break_score_vector[i] += break_score;

                // normalised break score by break counts
                double norm_break_count = (double)break_count/(double)total_breaks;
                double norm_break_score = prob*norm_break_count;
                norm_break_score_vector[i] += norm_break_score;

                // reset kmer break counter
                it->second = std::make_pair(std::get<0>(it->second), 0);       
            }
        }
        nr_of_breaks_vector[i] = total_breaks;
        
        // normalised break score by sequence length
        double norm_by_len = (double)break_score_vector[i]/(double)path[i].size();
        break_score_norm_by_len_vector[i] = norm_by_len;
    }

    // sort the break_score_vector in descending order
    // initialise original index locations
    std::vector<size_t> idx(break_score_vector.size());
    std::iota(idx.begin(), idx.end(), 0); 

    // sort indices based on comparing values in break_score_vector
    std::sort(idx.begin(), idx.end(), 
         [&break_score_vector](size_t idx_1, size_t idx_2){ 
            return break_score_vector[idx_1] > break_score_vector[idx_2]; 
    });

    // reorder elements
    std::vector<double> sorted_break_score(idx.size());
    std::vector<double> sorted_bp_score_norm_by_break_freqs(idx.size());
    std::vector<double> sorted_bp_score_norm_by_len(idx.size());
    std::vector<std::string> sorted_path(idx.size());
    std::vector<int> path_len(idx.size());
    std::vector<int> sorted_nr_of_breaks_vector(idx.size());
    std::vector<int> lev_dist_vs_true(idx.size());

    for(int i = 0; i < idx.size(); i++){
        sorted_break_score[i] = break_score_vector[idx[i]];
        sorted_bp_score_norm_by_break_freqs[i] = norm_break_score_vector[idx[i]];
        sorted_bp_score_norm_by_len[i] = break_score_norm_by_len_vector[idx[i]];
        sorted_path[i] = path[idx[i]];
        path_len[i] = path[idx[i]].length();
        sorted_nr_of_breaks_vector[i] = nr_of_breaks_vector[idx[i]];

        // calculate levenshtein distance of assembled vs true solution
        int lev_dist = calc_levenshtein(path[idx[i]], true_solution);
        lev_dist_vs_true[i] = lev_dist;
    }

    return Rcpp::List::create(
        Rcpp::Named("sequence") = sorted_path,
        Rcpp::Named("sequence_len") = path_len,
        Rcpp::Named("bp_score") = sorted_break_score,
        Rcpp::Named("bp_score_norm_by_break_freqs") = sorted_bp_score_norm_by_break_freqs,
        Rcpp::Named("bp_score_norm_by_len") = sorted_bp_score_norm_by_len,
        Rcpp::Named("kmer_breaks") = sorted_nr_of_breaks_vector,
        Rcpp::Named("lev_dist_vs_true") = lev_dist_vs_true
    );
}