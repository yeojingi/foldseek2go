#ifndef MMSEQS_BLAST2GO_H
#define MMSEQS_BLAST2GO_H

#include <map>
#include <vector>
#include <utility>

#include "GeneOntology.h"

class IndexReader;

namespace Blast2Go {

// Evidence-code-free BLAST2GO-style annotation rule: similarity term (BLAST2GO's
// "max.sim", supplied here as backtrace-derived positives %) + GO-hierarchy
// abstraction term (goWeight * merged children), filtered by cutoff, deduplicated
// by the true path rule (drop an ancestor if a more specific descendant also survives).
std::vector<std::pair<GoID, float>> annotate(
    const std::map<size_t, float>& simScoresForQuery, // target(tGoId) -> similarity(0~100)
    IndexReader* tGoDbr,
    const GeneOntology& go,
    float goWeight, float cutoff);

}

#endif
