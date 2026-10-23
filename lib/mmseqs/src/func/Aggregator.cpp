#include "Aggregator.h"
#include <map>
#include <cstring>
#include <algorithm>
#include "Parameters.h"
#include "IndexReader.h"
#include "Scorer.h"
#include "Util.h"
#include "GeneOntology.h"
#include "Blast2Go.h"

// _func 엔트리에서 GO term 목록 파싱 (newline 구분, 각 줄 첫 토큰만)
static void parseGoTerms(const char* data, size_t len, std::vector<std::string>& terms) {
    const char* ptr = data;
    const char* end = data + len;
    while (ptr < end) {
        const char* nl = (const char*)memchr(ptr, '\n', end - ptr);
        size_t lineLen = nl ? (size_t)(nl - ptr) : (size_t)(end - ptr);
        if (lineLen > 0) {
            const char* sp = (const char*)memchr(ptr, ' ', lineLen);
            size_t goLen = sp ? (size_t)(sp - ptr) : lineLen;
            terms.push_back(std::string(ptr, goLen));
        }
        if (!nl) break;
        ptr = nl + 1;
    }
}

std::vector<std::pair<std::string, float>> Aggregator::aggregateOneQuery(
    const std::map<size_t, float>* scores,
    const std::map<size_t, float>* simScores,
    IndexReader* tGoDbr,
    const GeneOntology* go,
    unsigned int thread_idx)
{
    std::vector<std::pair<std::string, float>> result;
    int policy = Parameters::getInstance().policy;

    if (policy == 2) {
        if (simScores == NULL || go == NULL) {
            return result;
        }
        Parameters& par = Parameters::getInstance();
        std::vector<std::pair<GoID, float>> annotated =
            Blast2Go::annotate(*simScores, tGoDbr, *go, par.goWeight, par.annotationCutoff);
        for (size_t i = 0; i < annotated.size(); i++) {
            char buf[16];
            snprintf(buf, sizeof(buf), "GO:%07d", annotated[i].first);
            result.push_back(std::make_pair(std::string(buf), annotated[i].second));
        }
    } else if (policy == 0) {
        // best e-value target의 모든 GO → score 1.0
        size_t bestTarget = UINT_MAX;
        float minEvalue = std::numeric_limits<float>::max();
        for (std::map<size_t, float>::const_iterator it = scores->begin(); it != scores->end(); ++it) {
            if (it->second < minEvalue) {
                minEvalue = it->second;
                bestTarget = it->first;
            }
        }
        if (bestTarget == UINT_MAX) return result;

        char* goData = tGoDbr->sequenceReader->getData(bestTarget, thread_idx);
        size_t goLen  = tGoDbr->sequenceReader->getSeqLen(bestTarget);

        std::vector<std::string> terms;
        parseGoTerms(goData, goLen, terms);
        for (size_t i = 0; i < terms.size(); i++) {
            result.push_back(std::make_pair(terms[i], 1.0f));
        }

    } else if (policy == 3) {
        // k-NN voting: take the k nearest neighbors (lowest e-value), vote among just those
        int k = Parameters::getInstance().knnK;
        std::vector<std::pair<float, size_t>> byEvalue;
        for (std::map<size_t, float>::const_iterator it = scores->begin(); it != scores->end(); ++it) {
            byEvalue.push_back(std::make_pair(it->second, it->first));
        }
        std::sort(byEvalue.begin(), byEvalue.end());
        size_t kUsed = std::min((size_t)k, byEvalue.size());
        if (kUsed == 0) return result;

        std::map<std::string, unsigned int> counts;
        for (size_t i = 0; i < kUsed; i++) {
            size_t targetIdx = byEvalue[i].second;
            char* goData = tGoDbr->sequenceReader->getData(targetIdx, thread_idx);
            size_t goLen  = tGoDbr->sequenceReader->getSeqLen(targetIdx);
            std::vector<std::string> terms;
            parseGoTerms(goData, goLen, terms);
            for (size_t j = 0; j < terms.size(); j++) {
                counts[terms[j]]++;
            }
        }
        for (std::map<std::string, unsigned int>::const_iterator it = counts.begin(); it != counts.end(); ++it) {
            float score = (float)it->second / (float)kUsed;
            result.push_back(std::make_pair(it->first, score));
        }

    } else if (policy == 4) {
        // DiamondScore (DeepGOPlus / NetGO BLAST-KNN):
        // S(q,t) = Σ_h bitscore(q,h)·1[t∈GO(h)] / Σ_h bitscore(q,h)
        // simScores carries bitscore per hit (populated by ScoreAlignmentsBitscore)
        if (simScores == nullptr) return result;
        double totalBs = 0.0;
        for (auto& kv : *simScores) totalBs += kv.second;
        if (totalBs == 0.0) return result;

        std::map<std::string, double> goWeightSum;
        for (auto& kv : *simScores) {
            size_t tIdx = kv.first;
            double bs   = kv.second;
            char*  goData = tGoDbr->sequenceReader->getData(tIdx, thread_idx);
            size_t goLen  = tGoDbr->sequenceReader->getSeqLen(tIdx);
            std::vector<std::string> terms;
            parseGoTerms(goData, goLen, terms);
            // deduplicate within this hit so each GO counted once per hit
            std::sort(terms.begin(), terms.end());
            terms.erase(std::unique(terms.begin(), terms.end()), terms.end());
            for (auto& t : terms) goWeightSum[t] += bs;
        }
        for (auto& kv : goWeightSum) {
            result.push_back({kv.first, (float)(kv.second / totalBs)});
        }

    } else {
        // voting: 모든 target의 GO term 카운팅 → score = count / total
        std::map<std::string, unsigned int> counts;
        unsigned int total = 0;
        for (std::map<size_t, float>::const_iterator it = scores->begin(); it != scores->end(); ++it) {
            char* goData = tGoDbr->sequenceReader->getData(it->first, thread_idx);
            size_t goLen  = tGoDbr->sequenceReader->getSeqLen(it->first);
            std::vector<std::string> terms;
            parseGoTerms(goData, goLen, terms);
            for (size_t i = 0; i < terms.size(); i++) {
                counts[terms[i]]++;
            }
            total++;
        }
        if (total == 0) return result;
        for (std::map<std::string, unsigned int>::const_iterator it = counts.begin(); it != counts.end(); ++it) {
            float score = (float)it->second / (float)total;
            result.push_back(std::make_pair(it->first, score));
        }
    }

    return result;
}

void Aggregator::aggregateAll(
    const EvidenceScore& evidenceScores,
    const EvidenceScore* simScores,
    const IndexReader* qDbrHeader,
    const IndexReader* tGoDbr,
    const GeneOntology* go,
    unsigned int thread_idx,
    FILE* resultFP,
    FILE* formattedIdsFP,
    int formatMode)
{
    if (formatMode == 0) {
        fprintf(resultFP, "query\tprediction\tscore\n");
    }

    typedef std::map<unsigned int, std::map<size_t, float>>::const_iterator EvidenceIt;
    bool firstRow = true;
    for (EvidenceIt qit = evidenceScores.begin(); qit != evidenceScores.end(); ++qit) {
        unsigned int queryId = qit->first;
        const std::map<size_t, float>& scores = qit->second;
        const std::map<size_t, float>* querySimScores = (simScores != NULL) ? simScores->getQueryScores(queryId) : NULL;

        size_t qHeaderId = qDbrHeader->sequenceReader->getId(queryId);
        const char* qHeader = qDbrHeader->sequenceReader->getData(qHeaderId, thread_idx);
        std::string queryIdStr = Util::parseFastaHeader(qHeader);

        std::vector<std::pair<std::string, float>> goScores =
            aggregateOneQuery(&scores, querySimScores, (IndexReader*)tGoDbr, go, thread_idx);

        // Mode 0 and 1: write to resultFP (TSV or JSON)
        if (formatMode == 1) {
            // JSON row: {"query":"...","go":"...","score":0.00}
            for (size_t i = 0; i < goScores.size(); i++) {
                if (!firstRow) fprintf(resultFP, ",\n");
                fprintf(resultFP, "{\"query\":\"%s\",\"go\":\"%s\",\"score\":%.3f}",
                    queryIdStr.c_str(),
                    goScores[i].first.c_str(),
                    goScores[i].second);
                firstRow = false;
            }
        } else if (formatMode == 0) {
            // TSV format (mode 0)
            if (goScores.empty()) {
                fprintf(resultFP, "%s\t\n", queryIdStr.c_str());
                continue;
            }
            for (size_t i = 0; i < goScores.size(); i++) {
                fprintf(resultFP, "%s\t%s\t%.2f\n",
                    queryIdStr.c_str(),
                    goScores[i].first.c_str(),
                    goScores[i].second);
            }
        }
        
        // Mode 2: write formatted_ids output
        // query_id<tab>target_ids_semicolon_separated<tab>scores_semicolon_separated
        if ((formatMode == 2 || formattedIdsFP) && !scores.empty()) {
            FILE* targetFP = (formatMode == 2) ? resultFP : formattedIdsFP;
            if (targetFP) {
                fprintf(targetFP, "%s", queryIdStr.c_str());
                
                // Write all target IDs separated by semicolon for field2
                fprintf(targetFP, "\t");
                int targetIdx = 0;
                for (std::map<size_t, float>::const_iterator it = scores.begin(); it != scores.end(); ++it) {
                    if (targetIdx > 0) fprintf(targetFP, ";");
                    fprintf(targetFP, "%zu", it->first);
                    targetIdx++;
                }
                
                // Write all scores separated by semicolon for field3
                fprintf(targetFP, "\t");
                int scoreIdx = 0;
                for (std::map<size_t, float>::const_iterator it = scores.begin(); it != scores.end(); ++it) {
                    if (scoreIdx > 0) fprintf(targetFP, ";");
                    fprintf(targetFP, "%.2f", it->second);
                    scoreIdx++;
                }
                
                fprintf(targetFP, "\n");
            }
        }
    }
}
