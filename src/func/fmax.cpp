#include "Parameters.h"
#include "IndexReader.h"
#include "GeneOntology.h"
#include "Util.h"
#include "Debug.h"

#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <string>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <climits>
#include <algorithm>
#include <functional>

// ── helpers ──────────────────────────────────────────────────────────────────

static inline int goStrToId(const char* s, size_t len) {
    if (len > 3 && s[0] == 'G' && s[1] == 'O' && s[2] == ':')
        return std::atoi(s + 3);
    return -1;
}

// BFS up the GO DAG; returns propagated set (seeds included)
static std::unordered_set<int> propagate(const std::unordered_set<int>& seeds,
                                          const GeneOntology& go) {
    std::unordered_set<int> result = seeds;
    std::queue<int> q;
    for (int g : seeds) q.push(g);
    while (!q.empty()) {
        int cur = q.front(); q.pop();
        const GoNode* node = go.getGo(cur);
        if (!node) continue;
        for (int p : node->parentGoIds)
            if (result.insert(p).second) q.push(p);
    }
    return result;
}

static size_t intersectSize(const std::unordered_set<int>& a,
                             const std::unordered_set<int>& b) {
    size_t cnt = 0;
    for (int g : a) cnt += b.count(g);
    return cnt;
}

// ── per-namespace GoSet stored per query ─────────────────────────────────────

struct GoSet {
    std::unordered_set<int> bp, mf, cc;
    bool empty() const { return bp.empty() && mf.empty() && cc.empty(); }
    const std::unordered_set<int>& ns(const std::string& n) const {
        if (n == "biological_process") return bp;
        if (n == "molecular_function") return mf;
        return cc;
    }
};

// ── F-max for one namespace ───────────────────────────────────────────────────

struct FmaxResult { float fmax, threshold, prec, rec; };

static FmaxResult computeFmax(
    const std::string& nsName,
    const std::unordered_map<std::string, GoSet>& groundTruth,
    const std::unordered_map<std::string, std::unordered_map<int,float>>& preds,
    const GeneOntology& go)
{
    // n = queries with at least one GT term in this namespace
    int n = 0;
    for (auto& [e, gs] : groundTruth)
        if (!gs.ns(nsName).empty()) n++;
    if (n == 0) return {0,0,0,0};

    FmaxResult best{0, 0, 0, 0};

    for (int ti = 0; ti <= 100; ti++) {
        float t = ti / 100.0f;
        double sumP = 0; int nP = 0;
        double sumR = 0;

        for (auto& [entry, gs] : groundTruth) {
            const auto& G = gs.ns(nsName);
            if (G.empty()) continue;

            auto pit = preds.find(entry);
            if (pit == preds.end()) { sumR += 0.0; continue; }

            // seeds: predicted terms in this namespace with score >= t
            std::unordered_set<int> seeds;
            for (auto& [gid, score] : pit->second) {
                if (score < t) continue;
                const GoNode* node = go.getGo(gid);
                if (node && node->category == nsName) seeds.insert(gid);
            }
            if (seeds.empty()) { sumR += 0.0; continue; }

            std::unordered_set<int> P = propagate(seeds, go);
            size_t inter = intersectSize(G, P);

            sumP += (double)inter / P.size();
            nP++;
            sumR += (double)inter / G.size();
        }

        if (nP == 0) continue;
        double prec = sumP / nP;
        double rec  = sumR / n;
        if (prec + rec == 0.0) continue;
        float F = (float)(2.0 * prec * rec / (prec + rec));
        if (F > best.fmax) { best.fmax = F; best.threshold = t;
                              best.prec = (float)prec; best.rec = (float)rec; }
    }
    return best;
}

// ── main entry point ──────────────────────────────────────────────────────────

int fmax(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    const int thread_idx = 0;
    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);

    // par.db1 = queryDB  (ground truth in _func, GOG in _func_gog)
    // par.db2 = prediction TSV  (query \t GO \t score)
    // par.db3 = output file

    GeneOntology go(par.db1 + "_func_gog");

    // ── load ground truth ────────────────────────────────────────────────────
    IndexReader qGoDbr(par.db1, par.threads, IndexReader::GO,
        touch ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0,
        DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    IndexReader qDbrHeader(par.db1, par.threads, IndexReader::SRC_HEADERS,
        touch ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);

    std::unordered_map<std::string, GoSet> groundTruth;
    groundTruth.reserve(qGoDbr.sequenceReader->getSize());

    for (size_t i = 0; i < qGoDbr.sequenceReader->getSize(); i++) {
        unsigned int key = qGoDbr.sequenceReader->getDbKey(i);
        char* data = qGoDbr.sequenceReader->getData(i, thread_idx);
        size_t len  = qGoDbr.sequenceReader->getSeqLen(i);

        size_t hid = qDbrHeader.sequenceReader->getId(key);
        if (hid == (size_t)UINT_MAX) continue;
        const char* hdr = qDbrHeader.sequenceReader->getData(hid, thread_idx);
        std::string entry = Util::parseFastaHeader(hdr);

        // parse GO terms
        std::unordered_set<int> seeds;
        const char* ptr = data, *end = data + len;
        while (ptr < end) {
            const char* nl = (const char*)memchr(ptr, '\n', end - ptr);
            size_t lineLen = nl ? (size_t)(nl - ptr) : (size_t)(end - ptr);
            if (lineLen > 0) {
                const char* sp = (const char*)memchr(ptr, ' ', lineLen);
                size_t goLen = sp ? (size_t)(sp - ptr) : lineLen;
                int gid = goStrToId(ptr, goLen);
                if (gid >= 0) seeds.insert(gid);
            }
            if (!nl) break;
            ptr = nl + 1;
        }
        if (seeds.empty()) continue;

        std::unordered_set<int> all = propagate(seeds, go);
        GoSet gs;
        for (int g : all) {
            const GoNode* n = go.getGo(g);
            if (!n) continue;
            if      (n->category == "biological_process") gs.bp.insert(g);
            else if (n->category == "molecular_function")  gs.mf.insert(g);
            else if (n->category == "cellular_component")  gs.cc.insert(g);
        }
        if (!gs.empty()) groundTruth[entry] = std::move(gs);
    }

    Debug(Debug::INFO) << "Loaded ground truth for " << groundTruth.size() << " queries\n";

    // ── load predictions ──────────────────────────────────────────────────────
    // format: query \t GO:XXXXXXX \t score
    std::unordered_map<std::string, std::unordered_map<int,float>> preds;
    {
        FILE* fp = fopen(par.db2.c_str(), "r");
        if (!fp) {
            Debug(Debug::ERROR) << "Cannot open prediction file: " << par.db2 << "\n";
            EXIT(EXIT_FAILURE);
        }
        char buf[4096];
        while (fgets(buf, sizeof(buf), fp)) {
            if (strncmp(buf, "query\t", 6) == 0) continue;  // header
            char* t1 = strchr(buf, '\t');
            if (!t1) continue;
            *t1 = '\0';
            char* t2 = strchr(t1+1, '\t');
            if (!t2) continue;
            *t2 = '\0';
            std::string entry(buf);
            int gid = goStrToId(t1+1, strlen(t1+1));
            if (gid < 0) continue;
            float score = (float)atof(t2+1);
            preds[entry][gid] = std::max(preds[entry][gid], score);
        }
        fclose(fp);
    }

    Debug(Debug::INFO) << "Loaded predictions for " << preds.size() << " queries\n";

    // ── compute F-max per namespace ───────────────────────────────────────────
    static const char* NS[] = { "biological_process", "molecular_function", "cellular_component" };
    static const char* LABELS[] = { "BP", "MF", "CC" };

    FILE* out = fopen(par.db3.c_str(), "w");
    if (!out) {
        Debug(Debug::ERROR) << "Cannot open output file: " << par.db3 << "\n";
        EXIT(EXIT_FAILURE);
    }
    fprintf(out, "namespace\tFmax\tthreshold\tprecision\trecall\n");

    for (int k = 0; k < 3; k++) {
        FmaxResult r = computeFmax(NS[k], groundTruth, preds, go);
        fprintf(out, "%s\t%.4f\t%.2f\t%.4f\t%.4f\n",
                LABELS[k], r.fmax, r.threshold, r.prec, r.rec);
        fprintf(stdout, "%s  Fmax=%.4f  t=%.2f  P=%.4f  R=%.4f\n",
                LABELS[k], r.fmax, r.threshold, r.prec, r.rec);
    }
    fclose(out);

    return EXIT_SUCCESS;
}
