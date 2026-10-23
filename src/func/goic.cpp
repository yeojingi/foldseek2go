#include "Parameters.h"
#include "IndexReader.h"
#include "GeneOntology.h"
#include "GoEnrichment.h"
#include "Debug.h"

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <climits>

// CAFA-style Information Accretion (Clark WT, Radivojac P: Information-theoretic
// evaluation of predicted ontological annotations. Bioinformatics, 29(13), i53-i61
// (2013), https://academic.oup.com/bioinformatics/article/29/13/i53/195366):
// IA(t) = -log2( P(t) / P(parents(t)) ), where
// P(parents(t)) is the probability that ALL of t's direct parents are annotated
// simultaneously (true path rule guarantees count(t) <= that joint count for every
// parent). N cancels out of the ratio, so this only needs raw counts:
//   IA(t) = -log2( count(t) / |intersection of postings(p) for p in parents(t)| )
// A term with no parents (an ontology root) gets IA = 0.

int goic(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    const unsigned int thread_idx = 0;
    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);

    GeneOntology go(par.db1 + "_func_gog");

    IndexReader* tGoDbr = new IndexReader(par.db1, par.threads, IndexReader::GO,
        touch ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0,
        DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);

    // _func(_gog) can be a much larger DB reused as-is (e.g. symlinked to a bigger
    // source DB's _func), so its index can carry entries that aren't actually part of
    // {db}. Only count entries that are also present in {db}.index.
    IndexReader baseDbr(par.db1, par.threads, IndexReader::SRC_SEQUENCES,
        touch ? IndexReader::PRELOAD_INDEX : 0,
        DBReader<unsigned int>::USE_INDEX);

    size_t N = tGoDbr->sequenceReader->getSize();

    // postings[t] = sorted sequence indices annotated with t (direct or propagated via
    // true path rule). Built in one pass over the background, indices ascending -> already
    // sorted, no separate sort step needed.
    std::unordered_map<int, std::vector<unsigned int>> postings;
    postings.reserve(60000);

    Debug(Debug::INFO) << "Building GO postings from background...\n";
    Debug::Progress progress(N);
    for (size_t i = 0; i < N; i++) {
        progress.updateProgress();
        unsigned int dbKey = tGoDbr->sequenceReader->getDbKey(i);
        if (baseDbr.sequenceReader->getId(dbKey) == UINT_MAX) continue;

        char* data = tGoDbr->sequenceReader->getData(i, thread_idx);
        size_t len = tGoDbr->sequenceReader->getSeqLen(i);

        std::vector<int> directIds;
        GoEnrichment::parseGoIds(data, len, directIds);
        if (directIds.empty()) continue;

        std::unordered_set<int> induced(directIds.begin(), directIds.end());
        for (int gid : directIds) {
            for (int anc : go.getAncestors(gid)) {
                induced.insert(anc);
            }
        }

        for (int gid : induced) {
            postings[gid].push_back((unsigned int)i);
        }
    }
    Debug(Debug::INFO) << "Built postings for " << postings.size() << " GO terms\n";

    static const std::vector<unsigned int> EMPTY_POSTINGS;
    auto getPostings = [&postings](int id) -> const std::vector<unsigned int>& {
        std::unordered_map<int, std::vector<unsigned int>>::const_iterator it = postings.find(id);
        return (it != postings.end()) ? it->second : EMPTY_POSTINGS;
    };

    FILE* out = fopen(par.db2.c_str(), "w");
    if (!out) {
        Debug(Debug::ERROR) << "Cannot open output file: " << par.db2 << "\n";
        EXIT(EXIT_FAILURE);
    }
    fprintf(out, "go_id\tgo_name\tcategory\tcount\tparent_count\tia\n");

    for (const auto& kv : postings) {
        int gid = kv.first;
        const std::vector<unsigned int>& countT = kv.second;

        const GoNode* node = go.getGo(gid);
        const char* goName = node ? node->goName.c_str() : "unknown";
        const char* goCat  = node ? node->category.c_str() : "unknown";

        size_t parentCount;
        if (node == nullptr || node->parentGoIds.empty()) {
            // ontology root: no parents to condition on, IA = 0 by convention
            parentCount = countT.size();
        } else {
            std::vector<unsigned int> denom = getPostings(node->parentGoIds[0]);
            for (size_t pi = 1; pi < node->parentGoIds.size(); pi++) {
                const std::vector<unsigned int>& next = getPostings(node->parentGoIds[pi]);
                std::vector<unsigned int> tmp;
                tmp.reserve(std::min(denom.size(), next.size()));
                std::set_intersection(denom.begin(), denom.end(), next.begin(), next.end(),
                                       std::back_inserter(tmp));
                denom = std::move(tmp);
            }
            parentCount = denom.size();
        }

        float ia = (parentCount > 0)
            ? -std::log2((float)countT.size() / (float)parentCount)
            : 0.0f;

        char goIdBuf[16];
        snprintf(goIdBuf, sizeof(goIdBuf), "GO:%07d", gid);
        fprintf(out, "%s\t%s\t%s\t%zu\t%zu\t%.4f\n",
                goIdBuf, goName, goCat, countT.size(), parentCount, ia);
    }

    fclose(out);
    delete tGoDbr;

    Debug(Debug::INFO) << "Done. IA written for " << postings.size() << " GO terms to " << par.db2 << "\n";
    return EXIT_SUCCESS;
}
