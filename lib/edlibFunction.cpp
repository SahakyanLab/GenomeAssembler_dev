#include <cstdio>
#include "edlib/edlib.h"
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int Levenshtein(const char* seq1, const char* seq2){
  static std::string work1; static std::string work2;
  work1 = seq1; work2 = seq2;
  
  const char *query  = work1.c_str();
  const char *target = work2.c_str();
  
  EdlibAlignResult result = edlibAlign(
    query, strlen(query), target, strlen(target), 
    edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
  
  if (result.status == EDLIB_STATUS_OK) {
    return result.editDistance;
  }
  return 0;
}