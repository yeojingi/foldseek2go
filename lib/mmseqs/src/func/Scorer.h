#ifndef MMSEQS_SCORER_H
#define MMSEQS_SCORER_H
#include "DBReader.h"
#include "IndexReader.h"
#include <map>

class BaseMatrix;

class EvidenceScore {
public:
  EvidenceScore() {}
  ~EvidenceScore() {}
  void addScore(unsigned int query_id, size_t target_id, float score) {
    aln_scores[query_id][target_id] = score;
  }

  float getScore(unsigned int query_id, size_t target_id) {
    return aln_scores[query_id][target_id];
  }

  const std::map<size_t, float>* getQueryScores(unsigned int query_id) const {
    auto it = aln_scores.find(query_id);
    return (it != aln_scores.end()) ? &(it->second) : nullptr;
  }


  auto begin() { return aln_scores.begin(); }
  auto end() { return aln_scores.end(); }
  auto begin() const { return aln_scores.begin(); }
  auto end() const { return aln_scores.end(); }

  void print() const {
    for (const auto& qpair : aln_scores) {
        unsigned int query_id = qpair.first;
        const auto& targetMap = qpair.second;

        std::cout << "Query " << query_id << ":\n";
        for (const auto& tpair : targetMap) {
            unsigned int target_id = tpair.first;
            float score = tpair.second;
            std::cout << "   Target " << target_id 
                      << " → Score: " << score 
                      << "\n";
        }
        std::cout << "\n";
    }
  }
private:
  std::map<unsigned int, std::map<size_t, float> > aln_scores;
};

void ScoreAlignments(DBReader<unsigned int>* , EvidenceScore* , IndexReader* );

// Bitscore per hit (res.score), used for policy 4 (DiamondScore/DeepGOPlus).
void ScoreAlignmentsBitscore(DBReader<unsigned int>* alnDbr, EvidenceScore* bsScore, IndexReader* tGoDbr);

// Positives-style (BLAST2GO "similarity") score per hit, computed from the alignment
// backtrace and a substitution matrix. Requires alnDbr to have been created with -a
// (backtrace); exits with an error otherwise.
void ScoreAlignmentsSimilarity(DBReader<unsigned int>* alnDbr, EvidenceScore* simScore, IndexReader* tGoDbr, IndexReader* qDbr, IndexReader* tDbr, BaseMatrix* subMat);

#endif