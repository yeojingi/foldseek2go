#ifndef MMSEQS_AGGREGATOR_H
#define MMSEQS_AGGREGATOR_H
#include <map>
#include <vector>
#include <string>
#include "IndexReader.h"
#include "Scorer.h"

class GeneOntology;

class Aggregator {
public:
  Aggregator() {}
  ~Aggregator() {}

  // simScores/go are only required for policy 2 (BLAST2GO-style); pass nullptr for
  // policy 0/1.
  std::vector<std::pair<std::string, float>> aggregateOneQuery(
      const std::map<size_t, float>* scores,
      const std::map<size_t, float>* simScores,
      IndexReader* tGoDbr,
      const GeneOntology* go,
      unsigned int thread_idx);

  void aggregateAll(
      const EvidenceScore& evidenceScores,
      const EvidenceScore* simScores,
      const IndexReader* qDbrHeader,
      const IndexReader* tGoDbr,
      const GeneOntology* go,
      unsigned int thread_idx,
      FILE*,
      FILE*,
      int formatMode
  );
};

#endif
