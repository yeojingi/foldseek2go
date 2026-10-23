#include "NcbiTaxonomy.h"
#include "Parameters.h"
#include "DBWriter.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include "FastSort.h"
#include "MappingReader.h"
#include "IndexReader.h"
#include "FuncReader.h"
#include "Matcher.h"
#include "goparser.h"
#include "GeneOntology.h"
#include "Aggregator.h"
#include "Scorer.h"
#include "SubstitutionMatrix.h"
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <cstring>
#include <cstdlib>
#include <climits>
#include <unistd.h>
#include <map>

#include "function.html.h"

#ifdef OPENMP
#include <omp.h>
#endif

int functionreport(int argc, const char **argv, const Command &command) {

unsigned int thread_idx = 0;
#ifdef OPENMP
thread_idx = omp_get_thread_num();
#endif

    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    const int formatMode = par.funcFormatMode;
    const bool devMode = par.devMode;
    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
    const bool sameDB = par.db1.compare(par.db2) == 0 ? true : false;

    GeneOntology go(par.db2 + "_func_gog");

    // policy 2 needs actual sequence data (to compute backtrace-based similarity);
    // other policies only ever look at the sequence index.
    const int dbaccessMode = (par.policy == 2)
        ? (DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA)
        : DBReader<unsigned int>::USE_INDEX;
    IndexReader qDbr(par.db1, par.threads, IndexReader::SRC_SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0, dbaccessMode);
    IndexReader qDbrHeader(par.db1, par.threads, IndexReader::SRC_HEADERS, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);

    IndexReader *tDbr;
    IndexReader *tDbrHeader;
    if (sameDB) {
        tDbr = &qDbr;
        tDbrHeader = &qDbrHeader;
    } else {
        tDbr = new IndexReader(par.db2, par.threads, IndexReader::SRC_SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0, dbaccessMode);
        tDbrHeader = new IndexReader(par.db2, par.threads, IndexReader::SRC_HEADERS, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);
    }

    DBReader<unsigned int> alnDbr(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    alnDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    IndexReader* tGoDbr = new IndexReader(par.db2, par.threads, IndexReader::GO, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    std::string db2NoIndexName = PrefilteringIndexReader::dbPathWithoutIndex(par.db2);
    FuncReader* funcMapping = new FuncReader(db2NoIndexName);

    // dev-mode: open query's _func DB to read query GO annotations
    IndexReader* qGoDbr = nullptr;
    if (devMode) {
        qGoDbr = new IndexReader(par.db1, par.threads, IndexReader::GO, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    }

    EvidenceScore evidenceScore;
    ScoreAlignments(&alnDbr, &evidenceScore, tGoDbr);
    std::cout << "Scored" << std::endl;

    // policy 2: similarity score; policy 4: bitscore (DiamondScore)
    EvidenceScore simScore;
    SubstitutionMatrix* subMat = nullptr;
    if (par.policy == 2) {
        subMat = new SubstitutionMatrix(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, 0.0);
        ScoreAlignmentsSimilarity(&alnDbr, &simScore, tGoDbr, &qDbr, tDbr, subMat);
    } else if (par.policy == 4) {
        ScoreAlignmentsBitscore(&alnDbr, &simScore, tGoDbr);
    }

    if (formatMode == 2) {
        // lookup: query entry(string) -> numericId
        std::map<std::string, unsigned int> queryLookup;

        typedef std::map<unsigned int, std::map<size_t, float>>::const_iterator EvidenceIt;
        for (EvidenceIt qit = evidenceScore.begin(); qit != evidenceScore.end(); ++qit) {
            unsigned int queryNumId = qit->first;
            size_t qHeaderId = qDbrHeader.sequenceReader->getId(queryNumId);
            const char* qHeader = qDbrHeader.sequenceReader->getData(qHeaderId, thread_idx);
            std::string queryStr = Util::parseFastaHeader(qHeader);
            queryLookup[queryStr] = queryNumId;
        }

        // --- Write .html file ---
        auto jsEscape = [](const std::string& s) -> std::string {
            std::string out;
            out.reserve(s.size());
            for (size_t i = 0; i < s.size(); i++) {
                if (s[i] == '\\' || s[i] == '"') out += '\\';
                out += s[i];
            }
            return out;
        };

        std::string gogFilePath = par.db2 + "_func_gog";

        std::string gogUrlPath = gogFilePath;
        std::string db4Base = par.db4;
        if (db4Base.size() >= 5 && db4Base.substr(db4Base.size() - 5) == ".html")
            db4Base = db4Base.substr(0, db4Base.size() - 5);
        std::string htmlPath = db4Base + ".html";
        FILE* htmlFP = FileUtil::openAndDelete(htmlPath.c_str(), "w");

        // Find /* INJECT_DATA */ placeholder in function_html template
        static const char INJECT_MARKER[] = "/* INJECT_DATA */";
        const size_t MARKER_LEN = sizeof(INJECT_MARKER) - 1;
        const char* html = (const char*)function_html;
        size_t html_len = (size_t)function_html_len;
        const char* inject_pos = nullptr;
        for (size_t i = 0; i + MARKER_LEN <= html_len; i++) {
            if (memcmp(html + i, INJECT_MARKER, MARKER_LEN) == 0) {
                inject_pos = html + i;
                break;
            }
        }

        // Write template up to placeholder
        if (inject_pos) {
            fwrite(html, 1, (size_t)(inject_pos - html), htmlFP);
        }

        // paths — all absolute
        // Note: argv[0] here is the first DB arg (Application.cpp dispatches argv+2).
        // The real binary path is stored in the MMSEQS env var by Application.cpp.
        char mmseqsAbsBuf[PATH_MAX];
        const char* mmseqsEnv = getenv("MMSEQS");
        std::string mmseqsPath = mmseqsEnv ? std::string(mmseqsEnv) : "";
        if (!mmseqsPath.empty() && realpath(mmseqsPath.c_str(), mmseqsAbsBuf) != nullptr)
            mmseqsPath = std::string(mmseqsAbsBuf);

        char cwdBuf[PATH_MAX];
        std::string cwd = (getcwd(cwdBuf, PATH_MAX) != nullptr) ? std::string(cwdBuf) : ".";

        auto toAbs = [&](const std::string& p) -> std::string {
            if (!p.empty() && p[0] == '/') return p;
            return cwd + "/" + p;
        };

        // Inject dynamic data (PATHS, GOG_PATH, TOTAL_COUNT, IDS_PATH, ENRICH_DB_PATH, LCA_PATH)
        std::string enrichDbAbs = par.enrichDb.empty() ? "" : toAbs(par.enrichDb);
        std::string lcaPath = db4Base + "_lca";
        fprintf(htmlFP,
            "const PATHS = {query:\"%s\", target:\"%s\", aln:\"%s\", mmseqs:\"%s\", cwd:\"%s\"};\n"
            "const GOG_PATH = \"%s\";\n"
            "const TOTAL_COUNT = %zu;\n"
            "const IDS_PATH = \"%s\";\n"
            "const ENRICH_DB_PATH = \"%s\";\n"
            "const LCA_PATH = \"%s\";\n"
            "const DEV_MODE = %s;\n",
            jsEscape(toAbs(par.db1)).c_str(),
            jsEscape(toAbs(par.db2)).c_str(),
            jsEscape(toAbs(par.db3)).c_str(),
            jsEscape(mmseqsPath).c_str(),
            jsEscape(cwd).c_str(),
            jsEscape(gogUrlPath).c_str(),
            queryLookup.size(),
            jsEscape(toAbs(db4Base + "_ids")).c_str(),
            jsEscape(enrichDbAbs).c_str(),
            jsEscape(toAbs(lcaPath)).c_str(),
            devMode ? "true" : "false"
        );

        // Write rest of template after placeholder
        if (inject_pos) {
            const char* after = inject_pos + MARKER_LEN;
            fwrite(after, 1, html_len - (size_t)(after - html), htmlFP);
        }

        fclose(htmlFP);

        // write _ids: entry\tnumeric_id\tprotein_name[\tquery_go_semicolon_separated]
        std::string idsPath = db4Base + "_ids";
        FILE* idsFP = FileUtil::openAndDelete(idsPath.c_str(), "w");
        for (auto it = queryLookup.begin(); it != queryLookup.end(); ++it) {
            size_t qHeaderId = qDbrHeader.sequenceReader->getId(it->second);
            const char* qHeader = qDbrHeader.sequenceReader->getData(qHeaderId, thread_idx);
            std::string fullHeader(qHeader);
            std::string protName;
            size_t sp = fullHeader.find(' ');
            if (sp != std::string::npos) {
                protName = fullHeader.substr(sp + 1);
                while (!protName.empty() && (protName.back() == '\n' || protName.back() == '\r' || protName.back() == ' '))
                    protName.pop_back();
            }
            if (devMode && qGoDbr) {
                size_t goIdx = qGoDbr->sequenceReader->getId(it->second);
                std::string queryGos;
                if (goIdx != UINT_MAX) {
                    char* goData = qGoDbr->sequenceReader->getData(goIdx, thread_idx);
                    size_t goLen  = qGoDbr->sequenceReader->getSeqLen(goIdx);
                    const char* ptr = goData;
                    const char* end = goData + goLen;
                    bool first = true;
                    while (ptr < end) {
                        const char* nl = (const char*)memchr(ptr, '\n', end - ptr);
                        size_t lineLen = nl ? (size_t)(nl - ptr) : (size_t)(end - ptr);
                        if (lineLen > 0) {
                            const char* spc = (const char*)memchr(ptr, ' ', lineLen);
                            size_t goTermLen = spc ? (size_t)(spc - ptr) : lineLen;
                            if (!first) queryGos += ';';
                            queryGos.append(ptr, goTermLen);
                            first = false;
                        }
                        if (!nl) break;
                        ptr = nl + 1;
                    }
                }
                fprintf(idsFP, "%s\t%u\t%s\t%s\n", it->first.c_str(), it->second, protName.c_str(), queryGos.c_str());
            } else {
                fprintf(idsFP, "%s\t%u\t%s\n", it->first.c_str(), it->second, protName.c_str());
            }
        }
        fclose(idsFP);

        // write _lca: entry\tnumId=lcaRank;numId=lcaRank;...
        {
            NcbiTaxonomy* tax = NcbiTaxonomy::openTaxonomy(par.db2);
            if (tax != nullptr) {
                MappingReader qMapper(par.db1);
                MappingReader tMapper(par.db2);
                FILE* lcaFP = FileUtil::openAndDelete(lcaPath.c_str(), "w");
                for (EvidenceIt qit = evidenceScore.begin(); qit != evidenceScore.end(); ++qit) {
                    unsigned int queryNumId = qit->first;
                    TaxID qTaxId = (TaxID)qMapper.lookup(queryNumId);
                    size_t qHid = qDbrHeader.sequenceReader->getId(queryNumId);
                    const char* qHdr = qDbrHeader.sequenceReader->getData(qHid, thread_idx);
                    std::string qEntry = Util::parseFastaHeader(qHdr);
                    fprintf(lcaFP, "%s\t", qEntry.c_str());
                    bool first = true;
                    for (auto tit = qit->second.begin(); tit != qit->second.end(); ++tit) {
                        size_t targetInternalIdx = tit->first;
                        // convert internal index → DB key (what mmseqs view outputs)
                        unsigned int targetDbKey = tGoDbr->sequenceReader->getDbKey(targetInternalIdx);
                        TaxID tTaxId = (TaxID)tMapper.lookup(targetDbKey);
                        const char* rank = "no rank";
                        if (qTaxId != 0 && tTaxId != 0) {
                            TaxID lcaTaxId = tax->LCA(qTaxId, tTaxId);
                            if (lcaTaxId != 0) {
                                const TaxonNode* node = tax->taxonNode(lcaTaxId, false);
                                if (node != nullptr)
                                    rank = tax->getString(node->rankIdx);
                            }
                        }
                        if (!first) fprintf(lcaFP, ";");
                        fprintf(lcaFP, "%u=%s", targetDbKey, rank);
                        first = false;
                    }
                    fprintf(lcaFP, "\n");
                }
                fclose(lcaFP);
                delete tax;
            }
        }

        std::cout << "Done\n"
                  << "  HTML:  " << htmlPath << "\n"
                  << "  IDs:   " << idsPath << "\n"
                  << "To view, run: python3 /path/to/m2g/server.py 8080\n"
                  << "Then open:   http://localhost:8080/" << htmlPath << "\n";
        delete tGoDbr;
        if (qGoDbr) delete qGoDbr;
        delete funcMapping;
        if (!sameDB) {
            delete tDbr;
            delete tDbrHeader;
        }
        if (subMat) delete subMat;
        alnDbr.close();
        return EXIT_SUCCESS;
    }

    // Mode 0 and 1
    FILE *resultFP = FileUtil::openAndDelete(par.db4.c_str(), "w");
    std::string formattedIdsPath = par.db4 + "_formatted_ids";
    FILE *formattedIdsFP = FileUtil::openAndDelete(formattedIdsPath.c_str(), "w");

    if (formatMode == 1) {
        fwrite(function_html, 1, function_html_len, resultFP);
    }

    Aggregator aggregator;
    aggregator.aggregateAll(
        evidenceScore,
        (par.policy == 2 || par.policy == 4) ? &simScore : nullptr,
        &qDbrHeader,
        tGoDbr,
        &go,
        thread_idx,
        resultFP,
        formattedIdsFP,
        formatMode
    );

    if (formatMode == 1) {
        fprintf(resultFP, "]);</script>");
    }

    std::cout << "Done" << std::endl;
    if (resultFP) fclose(resultFP);
    if (formattedIdsFP) fclose(formattedIdsFP);
    if (subMat) delete subMat;
    delete tGoDbr;
    delete funcMapping;
    if (!sameDB) {
        delete tDbr;
        delete tDbrHeader;
    }
    alnDbr.close();
    return EXIT_SUCCESS;
}
