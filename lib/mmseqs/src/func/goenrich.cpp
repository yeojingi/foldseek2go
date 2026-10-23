#include "Parameters.h"
#include "DBWriter.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include "IndexReader.h"
#include "Matcher.h"
#include "GeneOntology.h"
#include "Scorer.h"
#include "GoEnrichment.h"

#include <map>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <climits>

#ifdef OPENMP
#include <omp.h>
#endif

// ─── helpers ────────────────────────────────────────────────────────────────

static void buildGoIdStr(int id, char* buf, size_t bufLen) {
    snprintf(buf, bufLen, "GO:%07d", id);
}

// ─── command ────────────────────────────────────────────────────────────────

int goenrich(int argc, const char **argv, const Command &command) {
    unsigned int thread_idx = 0;
#ifdef OPENMP
    thread_idx = omp_get_thread_num();
#endif

    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
    (void)(par.db1 == par.db2); // sameDB unused for now

    GeneOntology go(par.db2 + "_func_gog");

    IndexReader qDbrHeader(par.db1, par.threads, IndexReader::SRC_HEADERS,
        touch ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);

    IndexReader *tGoDbr = new IndexReader(par.db2, par.threads, IndexReader::GO,
        touch ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0,
        DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);

    DBReader<unsigned int> alnDbr(par.db3.c_str(), par.db3Index.c_str(), par.threads,
        DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    alnDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    // ── build background GO counts ──────────────────────────────────────────
    Debug(Debug::INFO) << "Building background GO counts...\n";
    size_t bgTotal = tGoDbr->sequenceReader->getSize();
    std::unordered_map<int, unsigned int> bgGoCount; // goId -> # sequences in BG with it
    bgGoCount.reserve(60000);

    for (size_t i = 0; i < bgTotal; ++i) {
        char   *data = tGoDbr->sequenceReader->getData(i, thread_idx);
        size_t  len  = tGoDbr->sequenceReader->getSeqLen(i);
        std::vector<int> ids;
        GoEnrichment::parseGoIds(data, len, ids);
        // deduplicate per sequence
        std::sort(ids.begin(), ids.end());
        ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
        for (int gid : ids) bgGoCount[gid]++;
    }
    int N = (int)bgTotal;
    Debug(Debug::INFO) << "Background: " << N << " annotated sequences, "
                       << bgGoCount.size() << " distinct GO terms\n";

    // ── score alignments ────────────────────────────────────────────────────
    EvidenceScore evidenceScore;
    ScoreAlignments(&alnDbr, &evidenceScore, tGoDbr);

    // ── output: DBWriter indexed by query numeric ID ────────────────────────
    DBWriter enrichWriter(par.db4.c_str(), par.db4Index.c_str(),
                          1, par.compressed, Parameters::DBTYPE_GENERIC_DB);
    enrichWriter.open();

    // header as key UINT_MAX – skipped by viewer; write metadata to a sidecar instead
    std::string headerLine = "go_id\tgo_name\tcategory\tk\tM\tn\tN\tenrich_ratio\tpvalue\tadj_pvalue\n";

    typedef std::map<unsigned int, std::map<size_t, float>>::const_iterator EvidenceIt;
    size_t processed = 0;
    for (EvidenceIt qit = evidenceScore.begin(); qit != evidenceScore.end(); ++qit) {
        unsigned int queryId = qit->first;
        const std::map<size_t, float> &targets = qit->second;
        int M = (int)targets.size();
        if (M == 0) continue;

        // query string ID (for reference; not written to DB here)
        size_t qHdrIdx = qDbrHeader.sequenceReader->getId(queryId);
        const char *qHdr = qDbrHeader.sequenceReader->getData(qHdrIdx, thread_idx);
        std::string queryStr = Util::parseFastaHeader(qHdr);

        // count GO terms in aligned targets
        std::unordered_map<int, unsigned int> studyGoCount;
        studyGoCount.reserve(512);
        for (auto &tp : targets) {
            size_t tIdx = tp.first;
            char   *data = tGoDbr->sequenceReader->getData(tIdx, thread_idx);
            size_t  len  = tGoDbr->sequenceReader->getSeqLen(tIdx);
            std::vector<int> ids;
            GoEnrichment::parseGoIds(data, len, ids);
            std::sort(ids.begin(), ids.end());
            ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
            for (int gid : ids) studyGoCount[gid]++;
        }
        if (studyGoCount.empty()) continue;

        // compute p-values (enrichment: P(X>=k), depletion: P(X<=k))
        std::vector<int>    goList;
        std::vector<double> pvals;
        goList.reserve(studyGoCount.size());
        pvals.reserve(studyGoCount.size());
        for (auto &kv : studyGoCount) {
            int gid = kv.first;
            int k   = (int)kv.second;
            auto bgIt = bgGoCount.find(gid);
            int n = (bgIt != bgGoCount.end()) ? (int)bgIt->second : 0;
            goList.push_back(gid);
            double studyRate = (M > 0) ? (double)k / M : 0.0;
            double bgRate    = (N > 0 && n > 0) ? (double)n / N : 0.0;
            double pv = (studyRate >= bgRate)
                        ? GoEnrichment::hypergeomPval(N, n, M, k)
                        : GoEnrichment::hypergeomPvalLower(N, n, M, k);
            pvals.push_back(pv);
        }

        std::vector<double> adjPvals = GoEnrichment::bhCorrect(pvals);

        // sort by adj p-value ascending for the output
        std::vector<int> order(goList.size());
        for (int i = 0; i < (int)order.size(); ++i) order[i] = i;
        std::sort(order.begin(), order.end(),
                  [&](int a, int b){ return adjPvals[a] < adjPvals[b]; });

        // serialize entry
        std::string entry;
        entry.reserve(goList.size() * 80);
        char goIdBuf[16];
        char lineBuf[512];

        for (int oi : order) {
            int gid = goList[oi];
            int k   = (int)studyGoCount[gid];
            auto bgIt = bgGoCount.find(gid);
            int n = (bgIt != bgGoCount.end()) ? (int)bgIt->second : 0;
            double pval  = pvals[oi];
            double adj   = adjPvals[oi];
            double ratio = (n > 0 && M > 0)
                           ? ((double)k / M) / ((double)n / N)
                           : 0.0;

            const GoNode *node = go.getGo((GoID)gid);
            const char *goName = node ? node->goName.c_str()    : "unknown";
            const char *goCat  = node ? node->category.c_str()  : "unknown";
            buildGoIdStr(gid, goIdBuf, sizeof(goIdBuf));

            int len = snprintf(lineBuf, sizeof(lineBuf),
                "%s\t%s\t%s\t%d\t%d\t%d\t%d\t%.4f\t%.4e\t%.4e\n",
                goIdBuf, goName, goCat, k, M, n, N, ratio, pval, adj);
            if (len > 0) entry.append(lineBuf, len);
        }

        enrichWriter.writeData(entry.c_str(), entry.size(), queryId, 0);
        ++processed;
    }

    enrichWriter.close();

    // ── global enrichment: union of all study targets across ALL queries ────────
    Debug(Debug::INFO) << "Computing global enrichment (all queries combined)...\n";
    std::unordered_set<size_t> globalTargetSet;
    for (auto &qp : evidenceScore)
        for (auto &tp : qp.second)
            globalTargetSet.insert(tp.first);

    int globalM = (int)globalTargetSet.size();
    std::unordered_map<int, unsigned int> globalGoCount;
    globalGoCount.reserve(60000);
    for (size_t tIdx : globalTargetSet) {
        char   *data = tGoDbr->sequenceReader->getData(tIdx, thread_idx);
        size_t  len  = tGoDbr->sequenceReader->getSeqLen(tIdx);
        std::vector<int> ids;
        GoEnrichment::parseGoIds(data, len, ids);
        std::sort(ids.begin(), ids.end());
        ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
        for (int gid : ids) globalGoCount[gid]++;
    }

    if (!globalGoCount.empty()) {
        std::vector<int>    gGoList;
        std::vector<double> gPvals;
        gGoList.reserve(globalGoCount.size());
        gPvals.reserve(globalGoCount.size());
        for (auto &kv : globalGoCount) {
            int gid = kv.first, k = (int)kv.second;
            auto bgIt = bgGoCount.find(gid);
            int n = (bgIt != bgGoCount.end()) ? (int)bgIt->second : 0;
            gGoList.push_back(gid);
            double studyRate = (globalM > 0) ? (double)k / globalM : 0.0;
            double bgRate    = (N > 0 && n > 0) ? (double)n / N : 0.0;
            double pv = (studyRate >= bgRate)
                        ? GoEnrichment::hypergeomPval(N, n, globalM, k)
                        : GoEnrichment::hypergeomPvalLower(N, n, globalM, k);
            gPvals.push_back(pv);
        }
        std::vector<double> gAdj = GoEnrichment::bhCorrect(gPvals);
        std::vector<int> gOrder(gGoList.size());
        for (int i = 0; i < (int)gOrder.size(); ++i) gOrder[i] = i;
        std::sort(gOrder.begin(), gOrder.end(),
                  [&](int a, int b){ return gAdj[a] < gAdj[b]; });

        std::string globalPath = par.db4 + "_global";
        FILE *gfp = fopen(globalPath.c_str(), "w");
        if (gfp) {
            fprintf(gfp, "go_id\tgo_name\tcategory\tk\tM\tn\tN\tenrich_ratio\tpvalue\tadj_pvalue\n");
            char goIdBuf[16]; char lineBuf[512];
            for (int oi : gOrder) {
                int gid = gGoList[oi], k = (int)globalGoCount[gid];
                auto bgIt = bgGoCount.find(gid);
                int n = (bgIt != bgGoCount.end()) ? (int)bgIt->second : 0;
                double pval = gPvals[oi], adj = gAdj[oi];
                double ratio = (n > 0 && globalM > 0)
                               ? ((double)k / globalM) / ((double)n / N) : 0.0;
                const GoNode *node = go.getGo((GoID)gid);
                const char *goName = node ? node->goName.c_str()   : "unknown";
                const char *goCat  = node ? node->category.c_str() : "unknown";
                buildGoIdStr(gid, goIdBuf, sizeof(goIdBuf));
                int len = snprintf(lineBuf, sizeof(lineBuf),
                    "%s\t%s\t%s\t%d\t%d\t%d\t%d\t%.4f\t%.4e\t%.4e\n",
                    goIdBuf, goName, goCat, k, globalM, n, N, ratio, pval, adj);
                if (len > 0) fwrite(lineBuf, 1, len, gfp);
            }
            fclose(gfp);
            Debug(Debug::INFO) << "Global enrichment (" << gGoList.size() << " GO terms, "
                               << globalM << " targets) written to " << globalPath << "\n";
        }
    }

    delete tGoDbr;
    alnDbr.close();

    Debug(Debug::INFO) << "Done. Per-query enrichment for " << processed
                       << " queries written to " << par.db4 << "\n";
    Debug(Debug::INFO) << "  Columns: go_id, go_name, category, k, M, n, N, enrich_ratio, pvalue, adj_pvalue\n";
    return EXIT_SUCCESS;
}
