#ifndef MMSEQS_GENEONTOLOGY_H
#define MMSEQS_GENEONTOLOGY_H

// #include "StringBlock.h"

// #include <map>
#include <unordered_map>
#include <vector>
#include <string>

typedef int GoID;

struct GoNode {
public:
    int id;
    GoID goId;
    std::vector<GoID> parentGoIds;
    size_t nParents;
    std::string goName;
    std::string category; // biological_process, molecular_function, cellular_component
    // size_t rankIdx;
    // size_t nameIdx;

    GoNode() {};

    GoNode(int id, GoID goId, const std::vector<GoID>& parentGoIds, size_t nParents, const std::string &goName, const std::string &category = "")
    : id(id), goId(goId), parentGoIds(parentGoIds), nParents(nParents), goName(goName), category(category) {}
};


class GeneOntology {
public:
    static GeneOntology* openFuncDb(const std::string &database);
    GeneOntology(const std::string &fileName);
    ~GeneOntology();
    const GoNode* getGo(GoID id) const;
    const std::string getLineage(GoID id) const;
    // BFS of all ancestors (parents, grandparents, ...) of id, not including id itself
    std::vector<GoID> getAncestors(GoID id) const;

    // TaxonNode const * LCA(const std::vector<TaxID>& taxa) const;
    // TaxID LCA(TaxID taxonA, TaxID taxonB) const;
    // std::vector<std::string> AtRanks(TaxonNode const * node, const std::vector<std::string> &levels) const;
    // std::map<std::string, std::string> AllRanks(TaxonNode const *node) const;
    // std::string taxLineage(TaxonNode const *node, bool infoAsName = true);

    // static std::vector<std::string> parseRanks(const std::string& ranks);
    // static int findRankIndex(const std::string& rank);
    // static char findShortRank(const std::string& rank);

    // bool IsAncestor(TaxID ancestor, TaxID child);
    // TaxonNode const* taxonNode(TaxID taxonId, bool fail = true) const;
    // bool nodeExists(TaxID taxId) const;

    // std::unordered_map<TaxID, TaxonCounts> getCladeCounts(std::unordered_map<TaxID, unsigned int>& taxonCounts) const;

    // WeightedTaxResult weightedMajorityLCA(const std::vector<WeightedTaxHit> &setTaxa, const float majorityCutoff);

    // const char* getString(size_t blockIdx) const;

    // static std::pair<char*, size_t> serialize(const NcbiTaxonomy& taxonomy);
    // static NcbiTaxonomy* unserialize(char* data);

    // TaxonNode* taxonNodes;
    // size_t maxNodes;
private:
    // size_t loadNodes(std::vector<TaxonNode> &tmpNodes, const std::string &nodesFile);
    // size_t loadMerged(const std::string &mergedFile);
    // void loadNames(std::vector<TaxonNode> &tmpNodes, const std::string &namesFile);
    void loadGo(std::vector<GoNode> tmpNodes, const std::string &goFile);
    std::unordered_map<GoID, GoNode> goMap;
    // void elh(std::vector<std::vector<TaxID>> const & children, int node, int level, std::vector<int> &tmpE, std::vector<int> &tmpL);
    // void computeSparseTable();
    // int nodeId(TaxID taxId) const;

    // int RangeMinimumQuery(int i, int j) const;
    // int lcaHelper(int i, int j) const;

    // NcbiTaxonomy(TaxonNode* taxonNodes, size_t maxNodes, int maxTaxID, int *D, int *E, int *L, int *H, int **M, StringBlock<unsigned int> *block)
    //     : taxonNodes(taxonNodes), maxNodes(maxNodes), maxTaxID(maxTaxID), D(D), E(E), L(L), H(H), M(M), block(block), externalData(true), mmapData(NULL), mmapSize(0) {};
    // int maxTaxID;
    // int *D; // maps from taxID to node ID in taxonNodes
    // int *E; // for Euler tour sequence (size 2N-1)
    // int *L; // Level of nodes in tour sequence (size 2N-1)
    // int *H;
    // int **M;
    // StringBlock<unsigned int>* block;

    // bool externalData;
    // char* mmapData;
    // size_t mmapSize;

    // static const int SERIALIZATION_VERSION;
};

#endif
