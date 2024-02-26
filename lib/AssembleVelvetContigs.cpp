// [[Rcpp::plugins("cpp17")]]
#include <Rcpp.h>
#include <algorithm>
#include <numeric>
#include <random>

// fast hash map
#include <gtl/include/gtl/phmap.hpp>

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
 * Assemble each batch of contigs through brute force alignment.
 * @param contig_matrix list of shuffled contigs.
 * @param dbg_kmer size of kmer per de bruijn graph node.
 * @return top_5_percent_matrix top 5% of assembled solutions by length.
*/
// [[Rcpp::export]]
std::vector<std::string> assemble_contigs(
    const std::vector<std::string> &velvet_contigs, 
    const int &dbg_kmer, 
    const int &seed){
 
    // Randomly shuffle the order of contigs for assembly
    int matrix_size = 20000;
    std::mt19937 engine(seed);
    std::vector<std::vector<std::string>> contig_matrix;
    contig_matrix.resize(matrix_size);
    for(int i = 0; i < matrix_size; i++){
        std::vector<std::string> contigs_copy = velvet_contigs;
        std::shuffle(contigs_copy.begin(), contigs_copy.end(), engine);
        contig_matrix[i] = contigs_copy;
    }

    // scaffold matrix
    std::vector<std::vector<std::string>> scaffold_matrix;
    scaffold_matrix.resize(contig_matrix.size());

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
        scaffold_matrix[contig_ind] = contigs;

        // update total number of contigs
        contig_count += contigs.size();
    }

    // flatten list into character vector and discard duplicates
    std::vector<std::string> flat_contig_matrix(contig_count);
    int ind = 0;
    for(int i = 0; i < scaffold_matrix.size(); i++){
        std::vector<std::string> current_contigs = scaffold_matrix[i];
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

    return flat_contig_matrix;
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

    // hash map of a de novo sequence solution and vector of probabilites
    gtl::flat_hash_map<std::string, std::vector<double>> path_prob_dist;
    for(int i = 0; i < path.size(); i++){
        std::string current_path = path[i];

        // init vector of prob for current path
        std::vector<double> prob_dist(current_path.length()-kmer+1);

        // get probability values in rolling window
        for(int pos = 0; pos <= current_path.length()-kmer; pos++){
            std::string current_kmer = current_path.substr(pos, kmer);
            prob_dist[pos] = bp_matrix[current_kmer].first;
        }

        // push into hash map
        path_prob_dist[path[i]] = prob_dist;
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

    for(int i = 0; i < path.size(); i++){
        // init total count of kmer breaks
        int total_breaks = 0;

        for(const auto &kv : read_matrix){
            // extract key-value pair
            std::string read = kv.first; 
            int read_count = kv.second;

            // find exact match
            std::size_t pos = path[i].find(read);

            if(pos != std::string::npos){
                // get start pos of broken kmer
                // take max of start pos of kmer or start of path
                // to prevent having negative start positions for expanding into k-mers.
                int start_pos_ind = std::max(0, (int)pos-(kmer/2));
                int expanded_kmer = 8;

                // if start_pos_ind is zero, expand above into smaller k-mers.
                if(start_pos_ind == 0){
                    if(pos == 1){
                        // 2-mer
                        expanded_kmer = 2;
                    } else if(pos == 2){
                        // 4-mer
                        expanded_kmer = 4;
                    } else if(pos == 3){
                        // 6-mer
                        expanded_kmer = 6;
                    }
                }
                // don't need to check the same on the other end of the path
                // because we only take the breakpoint from the start of the read

                // extract broken kmer
                std::string broken_kmer = path[i].substr(start_pos_ind, expanded_kmer);

                // update counter in matrix
                bp_matrix[broken_kmer].second += read_count;
                total_breaks += read_count;
            }
        }

        // Loop over the bp_matrix
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
    std::vector<std::vector<double>> sorted_path_prob_dist(idx.size());

    for(int i = 0; i < idx.size(); i++){
        sorted_break_score[i] = break_score_vector[idx[i]];
        sorted_bp_score_norm_by_break_freqs[i] = norm_break_score_vector[idx[i]];
        sorted_bp_score_norm_by_len[i] = break_score_norm_by_len_vector[idx[i]];
        sorted_path[i] = path[idx[i]];
        sorted_path_len[i] = path_len[idx[i]];
        sorted_nr_of_breaks_vector[i] = nr_of_breaks_vector[idx[i]];
        sorted_path_prob_dist[i] = path_prob_dist[sorted_path[i]];

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
        Rcpp::Named("lev_dist_vs_true") = lev_dist_vs_true,
        Rcpp::Named("path_prob_dist") = Rcpp::wrap(sorted_path_prob_dist)
    );
}