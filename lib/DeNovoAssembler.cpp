// [[Rcpp::plugins("cpp17")]]
#include <Rcpp.h>
#include <algorithm>
#include <numeric>
#include <random>

// fast hash map
#include <gtl/include/gtl/phmap.hpp>

// parallel for loop
#include </usr/local/Cellar/libomp/15.0.7/include/omp.h>

// edlib dependency
#include "../lib/edlib/edlib.h"

using namespace Rcpp;

/**
 * Get reverse complement of a DNA string.
 * @param sequence kmer broken after aligning read with de novo sequence/
 * @return reverse_seq reverse complement of sequence.
*/
std::string reverse_complement(const std::string &sequence){
    // reverse complement map
    gtl::flat_hash_map<char, char> complement = {
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
int calc_levenshtein(const std::string &query, const std::string &target){
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
 * Remove duplicates
 * @param vec vector containing duplicates.
 * @return unique elements in vector.
*/
template <typename T>
void remove_duplicates(std::vector<T> &vec){
    std::sort(vec.begin(), vec.end());
    
    vec.erase(std::unique(
        vec.begin(), 
        vec.end()), 
        vec.end()
    );
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
gtl::flat_hash_map<int, std::vector<std::string>> get_contigs(
    const std::vector<std::string> &read_kmers, 
    const int &dbg_kmer, 
    const int &seed){    
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
    gtl::flat_hash_map<std::string, std::vector<std::string>> dict;

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
    gtl::flat_hash_map<std::string, std::pair<int, int>> balanced_count;

    // get unique prefix and suffix
    remove_duplicates(prefix);
    remove_duplicates(suffix);
    
    // combine results 
    prefix.insert(
        prefix.end(), 
        suffix.begin(), 
        suffix.end()
    );

    // init balanced counts to zeros
    for(int i = 0; i < prefix.size(); i++){
        std::string key = prefix[i];
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
            balanced_count[current_edge].first++;
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

    // remove duplicates
    remove_duplicates(contigs);

    // Randomly shuffle the order of contigs for assembly
    int matrix_size = 50000;
    std::mt19937 engine(seed);
    gtl::flat_hash_map<int, std::vector<std::string>> contig_matrix;
    contig_matrix.reserve(matrix_size);
    for(int i = 0; i < matrix_size; i++){
        std::vector<std::string> contigs_copy = contigs;
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
std::vector<std::string> assemble_contigs(
    const std::vector<std::string> &read_kmers, 
    const int &dbg_kmer, 
    const int &seed){

    // get contigs
    gtl::flat_hash_map<int, std::vector<std::string>> contig_matrix = get_contigs(
        read_kmers, dbg_kmer, seed
    );

    // init count of contigs
    int contig_count = 0;

    // Loop over all the contig subsets and return assemblies in-place.
    // Note. Below code looks worse than it appears. In reality, it executes fast
    for(int contig_ind = 0; contig_ind < contig_matrix.size(); contig_ind++){
        
        // get contigs
        std::vector<std::string> contigs = contig_matrix[contig_ind];

        for(int kmer = (dbg_kmer-1); kmer > 0; kmer--){
            bool len_changed = true;
            while(len_changed){
                int temp = contigs.size();
                for(int i = 0; i < contigs.size(); i++){
                    if (contigs[i] == "") continue;
                    for(int j = contigs.size()-1; j >= 0; j--){
                        if(contigs[i] != contigs[j]){
                            std::string suffix = contigs[i].substr(
                                contigs[i].size() - kmer, 
                                contigs[i].size()
                            );
                            std::string prefix = contigs[j].substr(0, kmer);

                            // If the suffix and prefix match, concatenate the two contigs
                            if(suffix == prefix){
                                std::string no_prefix = contigs[j].substr(
                                    kmer, 
                                    contigs[j].size()
                                );
                                contigs[i] = contigs[i].append(no_prefix);
                                contigs[j] = "";
                            }
                        }
                    }
                }
                for(int i = contigs.size()-1; i >= 0; i--){
                    if(contigs[i] == ""){
                        contigs.erase(contigs.begin()+i);
                    }
                }
                len_changed = temp != contigs.size();
            }
        }
        // Update the contig subset with the assembled contigs
        contig_matrix[contig_ind] = contigs;

        // update total number of contigs
        contig_count += contigs.size();
    }

    // flatten list into character vector and discard duplicates
    std::vector<std::string> flat_contig_matrix(contig_count);
    int ind = 0;
    for(const auto &kv : contig_matrix){
        std::vector<std::string> current_contigs = kv.second;
        for(int j = 0; j < current_contigs.size(); j++){
            flat_contig_matrix[ind] = current_contigs[j];
            ind++;
        }
    }

    // Discard duplicates
    remove_duplicates(flat_contig_matrix);
    
    // sort vector in descending order by length
    std::sort(
        flat_contig_matrix.begin(),
        flat_contig_matrix.end(),
        [](const std::string &seq_a, const std::string &seq_b) -> bool {
            return seq_a.length() > seq_b.length();
    });

    // Retain top 5% of assembled solution by length
    int top_5_percent = flat_contig_matrix.size()*0.05;
    std::vector<std::string> top_5_percent_matrix(
        flat_contig_matrix.begin(),
        flat_contig_matrix.begin()+top_5_percent
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
Rcpp::List calc_breakscore(const std::vector<std::string> &path, 
                           const std::vector<std::string> &sequencing_reads, 
                           const std::string &true_solution,
                           const int &kmer, 
                           const std::vector<std::string> &bp_kmer,
                           const std::vector<double> &bp_prob,
                           const int &num_threads){

    // hash map of kmers, probability values and break counts
    gtl::flat_hash_map<std::string, std::pair<double, int>> bp_matrix;
    for(int i = 0; i < bp_kmer.size(); i++){
        bp_matrix[bp_kmer[i]] = std::make_pair(bp_prob[i], 0);
    }

    // hash map of reads and number of counts
    gtl::flat_hash_map<std::string, int> read_matrix;
    for(const auto & read : sequencing_reads){
        read_matrix[read]++;
    }

    // init break score vector
    std::vector<double> break_score_vector(path.size());
    std::vector<double> norm_break_score_vector(path.size());
    std::vector<double> break_score_norm_by_len_vector(path.size());
    std::vector<double> nr_of_breaks_vector(path.size());
    std::vector<int> path_len(path.size());

    #pragma omp parallel for num_threads(num_threads) shared(read) private(i)
    {
        #pragma omp for 
        for(int i = 0; i < path.size(); i++){
            // init total count of kmer breaks
            int total_breaks = 0;

            #pragma omp critical
            for(const auto &kv : read_matrix){
                // extract key-value pair
                std::string read = kv.first; 
                int read_count = kv.second;

                // find exact match
                std::size_t pos = path[i].find(read);

                if(pos != std::string::npos){
                    // get start pos of broken kmer
                    // take max of start pos of kmer or start of path
                    int start_pos_ind = std::max(0, (int)pos-(kmer/2));

                    // extract broken kmer
                    std::string broken_kmer = path[i].substr(start_pos_ind, kmer);

                    // update counter in matrix
                    bp_matrix[broken_kmer].second += read_count;
                    total_breaks += read_count;
                }
            }

            // Loop over the bp_matrix
            #pragma omp for
            for(auto &kv : bp_matrix){
                double prob = kv.second.first;
                int break_count = kv.second.second;

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
                    kv.second.second = 0;
                }
            }
            nr_of_breaks_vector[i] = total_breaks;
            path_len[i] = (int)path[i].length();
            
            // normalised break score by sequence length
            double norm_by_len = (double)break_score_vector[i]/path_len[i];
            break_score_norm_by_len_vector[i] = norm_by_len;
        }
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
    std::vector<int> sorted_path_len(idx.size());
    std::vector<int> sorted_nr_of_breaks_vector(idx.size());
    std::vector<int> lev_dist_vs_true(idx.size());

    for(int i = 0; i < idx.size(); i++){
        sorted_break_score[i] = break_score_vector[idx[i]];
        sorted_bp_score_norm_by_break_freqs[i] = norm_break_score_vector[idx[i]];
        sorted_bp_score_norm_by_len[i] = break_score_norm_by_len_vector[idx[i]];
        sorted_path[i] = path[idx[i]];
        sorted_path_len[i] = path_len[idx[i]];
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