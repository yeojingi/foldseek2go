#include "Blast2Go.h"
#include "GoEnrichment.h"
#include "IndexReader.h"

#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <algorithm>

namespace Blast2Go {

std::vector<std::pair<GoID, float>> annotate(
    const std::map<size_t, float>& simScoresForQuery,
    IndexReader* tGoDbr,
    const GeneOntology& go,
    float goWeight, float cutoff)
{
    const unsigned int thread_idx = 0;

    // 1. DirectScore: best similarity per directly-annotated GO term
    std::unordered_map<GoID, float> directScore;
    for (const auto& kv : simScoresForQuery) {
        size_t targetIdx = kv.first;
        float similarity = kv.second;
        char* data = tGoDbr->sequenceReader->getData(targetIdx, thread_idx);
        size_t len = tGoDbr->sequenceReader->getSeqLen(targetIdx);
        std::vector<int> ids;
        GoEnrichment::parseGoIds(data, len, ids);
        for (int gid : ids) {
            auto it = directScore.find(gid);
            if (it == directScore.end() || similarity > it->second) {
                directScore[gid] = similarity;
            }
        }
    }
    if (directScore.empty()) {
        return std::vector<std::pair<GoID, float>>();
    }

    // 2. induced set = direct terms + all their ancestors (true path rule closure)
    std::unordered_set<GoID> induced;
    for (const auto& kv : directScore) {
        induced.insert(kv.first);
        for (GoID anc : go.getAncestors(kv.first)) {
            induced.insert(anc);
        }
    }

    // local reverse-parent map (only edges within the induced set)
    std::unordered_map<GoID, std::vector<GoID>> children;
    for (GoID n : induced) {
        const GoNode* node = go.getGo(n);
        if (!node) continue;
        for (GoID p : node->parentGoIds) {
            if (induced.count(p)) {
                children[p].push_back(n);
            }
        }
    }

    // 3. bottom-up inherited score: whatever supports a child also supports its parent
    std::unordered_map<GoID, float> inherited;
    std::function<float(GoID)> computeInherited = [&](GoID n) -> float {
        auto it = inherited.find(n);
        if (it != inherited.end()) {
            return it->second;
        }
        float best = 0.0f;
        auto dsIt = directScore.find(n);
        if (dsIt != directScore.end()) {
            best = dsIt->second;
        }
        auto chIt = children.find(n);
        if (chIt != children.end()) {
            for (GoID c : chIt->second) {
                float cs = computeInherited(c);
                if (cs > best) {
                    best = cs;
                }
            }
        }
        inherited[n] = best;
        return best;
    };
    for (GoID n : induced) {
        computeInherited(n);
    }

    // 4. abstraction term: reward nodes where multiple child branches converge
    std::unordered_map<GoID, float> finalScore;
    for (GoID n : induced) {
        auto chIt = children.find(n);
        float childCount = (chIt != children.end()) ? (float)chIt->second.size() : 0.0f;
        finalScore[n] = inherited[n] + goWeight * childCount;
    }

    // 5. cutoff
    std::unordered_set<GoID> survivors;
    for (const auto& kv : finalScore) {
        if (kv.second >= cutoff) {
            survivors.insert(kv.first);
        }
    }

    // 6. true path rule: drop a survivor if a more specific descendant also survived
    std::unordered_map<GoID, bool> hasSurvivorMemo;
    std::function<bool(GoID)> hasSurvivingDescendant = [&](GoID n) -> bool {
        auto it = hasSurvivorMemo.find(n);
        if (it != hasSurvivorMemo.end()) {
            return it->second;
        }
        bool result = false;
        auto chIt = children.find(n);
        if (chIt != children.end()) {
            for (GoID c : chIt->second) {
                if (survivors.count(c) || hasSurvivingDescendant(c)) {
                    result = true;
                    break;
                }
            }
        }
        hasSurvivorMemo[n] = result;
        return result;
    };

    std::vector<std::pair<GoID, float>> result;
    for (GoID n : survivors) {
        if (!hasSurvivingDescendant(n)) {
            result.push_back(std::make_pair(n, finalScore[n]));
        }
    }

    std::sort(result.begin(), result.end(),
              [](const std::pair<GoID, float>& a, const std::pair<GoID, float>& b) {
                  return a.second > b.second;
              });

    return result;
}

}
