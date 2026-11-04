#include "DBReader.h"
#include "IndexReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "LocalParameters.h"
#include "Matcher.h"
#include "Alignment.h"
#include "structureto3diseqdist.h"
#include "StructureSmithWaterman.h"
#include "StructureUtil.h"
#include "TMaligner.h"
#include "Coordinate16.h"
#include "LDDT.h"
#include "NcbiTaxonomy.h"
#include "MappingReader.h"
#include "PrefilteringIndexReader.h"

#include <cfloat>
#include <cstdint>

#ifdef OPENMP
#include <omp.h>
#endif

// need for sorting the results
static bool compareHitsByStructureBits(const Matcher::result_t &first, const Matcher::result_t &second) {
    if (first.score != second.score) {
        return first.score > second.score;
    }
    if (first.dbLen != second.dbLen) {
        return first.dbLen < second.dbLen;
    }
    return first.dbKey < second.dbKey;
}

// align2b exclusion: index of the first standard rank at or above the LCA of a and b,
// -1 if there is no standard rank below the root
static int lcaRankIndex(const NcbiTaxonomy *taxonomy, TaxID a, TaxID b) {
    TaxonNode const *node = taxonomy->taxonNode(taxonomy->LCA(a, b), false);
    while (node != NULL) {
        const std::string rank = taxonomy->getString(node->rankIdx);
        // newer taxdumps renamed "superkingdom" to "domain", and the Viruses node carries "acellular root"
        const int rankIdx = (rank == "domain" || rank == "acellular root")
                            ? NcbiTaxonomy::findRankIndex("superkingdom") : NcbiTaxonomy::findRankIndex(rank);
        if (rankIdx >= 0) {
            return rankIdx;
        }
        if (node->parentTaxId == node->taxId) {
            break;
        }
        node = taxonomy->taxonNode(node->parentTaxId, false);
    }
    return -1;
}

static bool isTooClose(int lcaRankIdx, int thresholdRankIdx) {
    return lcaRankIdx >= 0 && lcaRankIdx <= thresholdRankIdx;
}


static void structureAlignDefault(LocalParameters & par) {
    par.compBiasCorrectionScale = 0.5;
    par.alignmentType = LocalParameters::ALIGNMENT_TYPE_3DI_AA;
}

int alignStructure(StructureSmithWaterman & structureSmithWaterman,
                   StructureSmithWaterman & reverseStructureSmithWaterman,
                   Sequence & tSeqAA, Sequence & tSeq3Di,
                   unsigned int querySeqLen, unsigned int targetSeqLen,
                   EvalueNeuralNet & evaluer, std::pair<double, double> muLambda,
                   Matcher::result_t & res, std::string & backtrace,
                   Parameters & par, double evalThr, float covThr, int32_t minScore) {

    float seqId = 0.0;
    backtrace.clear();
    // align only score and end pos
    StructureSmithWaterman::s_align align = structureSmithWaterman.alignScoreEndPos<StructureSmithWaterman::PROFILE>(tSeqAA.numSequence, tSeq3Di.numSequence, targetSeqLen, par.gapOpen.values.aminoacid(),
                                                                                    par.gapExtend.values.aminoacid(), querySeqLen / 2);
    bool hasLowerCoverage = !(Util::hasCoverage(covThr, par.covMode, align.qCov, align.tCov));
    if(hasLowerCoverage){
        return -1;
    }
    // the reverse score is only subtracted, so the forward score is an upper bound
    if (static_cast<int32_t>(align.score1) < minScore) {
        return -1;
    }
    // we can already stop if this e-value isn't good enough, it wont be any better in the next step
    align.evalue = evaluer.computeEvalueCorr(align.score1, muLambda.first, muLambda.second);
    bool hasLowerEvalue = align.evalue > evalThr;
    if(hasLowerEvalue){
        return -1;
    }

    StructureSmithWaterman::s_align revAlign;
    //if(structureSmithWaterman.isProfileSearch()){
    //    revAlign.score1 = 0;
    //} else {
    revAlign = reverseStructureSmithWaterman.alignScoreEndPos<StructureSmithWaterman::PROFILE>(tSeqAA.numSequence, tSeq3Di.numSequence,
                                                                  targetSeqLen, par.gapOpen.values.aminoacid(),
                                                                  par.gapExtend.values.aminoacid(), querySeqLen / 2);
    //}
    int32_t score = static_cast<int32_t>(align.score1) - static_cast<int32_t>(revAlign.score1);
    if (score < minScore) {
        return -1;
    }
    align.evalue = evaluer.computeEvalueCorr(score, muLambda.first, muLambda.second);
    hasLowerEvalue = align.evalue > evalThr;
    if (hasLowerEvalue) {
        return -1;
    }

    bool blockAlignFailed = false;
    if (structureSmithWaterman.isProfileSearch() == false) {
        StructureSmithWaterman::s_align alignTmp = structureSmithWaterman.alignStartPosBacktraceBlock(
            tSeqAA.numSequence, tSeq3Di.numSequence, targetSeqLen, par.gapOpen.values.aminoacid(),
            par.gapExtend.values.aminoacid(), backtrace, align
        );

        if (align.score1 == UINT32_MAX) {
            Debug(Debug::WARNING) << "block-align failed, falling back to normal alignment\n";
            blockAlignFailed = true;
        } else {
            align = alignTmp;
        }
    }

    if (blockAlignFailed || structureSmithWaterman.isProfileSearch()) {
        align = structureSmithWaterman.alignStartPosBacktrace<StructureSmithWaterman::PROFILE>(tSeqAA.numSequence,
                                                                                               tSeq3Di.numSequence,
                                                                                               targetSeqLen,
                                                                                               par.gapOpen.values.aminoacid(),
                                                                                               par.gapExtend.values.aminoacid(),
                                                                                               par.alignmentMode,
                                                                                               backtrace, align,
                                                                                               par.covMode, covThr,
                                                                                               querySeqLen / 2);
    }

    unsigned int alnLength = Matcher::computeAlnLength(align.qStartPos1, align.qEndPos1, align.dbStartPos1, align.dbEndPos1);
    if(backtrace.size() > 0){
        alnLength = backtrace.size();
        seqId = Util::computeSeqId(par.seqIdMode, align.identicalAACnt, querySeqLen, targetSeqLen, alnLength);
    }
    align.score1 = score;
    res = Matcher::result_t(tSeqAA.getDbKey(), align.score1, align.qCov, align.tCov, seqId, align.evalue, alnLength,
                            align.qStartPos1, align.qEndPos1, querySeqLen, align.dbStartPos1, align.dbEndPos1, targetSeqLen, backtrace);
    return 0;
}


int computeAlternativeAlignment(StructureSmithWaterman & structureSmithWaterman,
                                StructureSmithWaterman & reverseStructureSmithWaterman,
                                Sequence & tSeqAA, Sequence & tSeq3Di,
                                unsigned int querySeqLen, unsigned int targetSeqLen,
                                EvalueNeuralNet & evaluer, std::pair<double, double> muLambda,
                                Matcher::result_t & result, Matcher::result_t & altRes,
                                std::string & backtrace, Parameters & par) {
    const unsigned char xAAIndex = tSeqAA.subMat->aa2num[static_cast<int>('X')];
    const unsigned char x3DiIndex = tSeq3Di.subMat->aa2num[static_cast<int>('X')];
    for (int pos = result.dbStartPos; pos < result.dbEndPos; ++pos) {
        tSeqAA.numSequence[pos] = xAAIndex;
        tSeq3Di.numSequence[pos] = x3DiIndex;
    }
    if (alignStructure(structureSmithWaterman, reverseStructureSmithWaterman,
                       tSeqAA, tSeq3Di, querySeqLen, targetSeqLen,
                       evaluer, muLambda, altRes, backtrace, par,
                       par.evalThr, par.covThr, INT32_MIN) == -1) {
        return -1;
    }
    if (Alignment::checkCriteria(altRes, false, par.evalThr, par.seqIdThr, par.alnLenThr, par.covMode, par.covThr)) {
        return 0;
    } else {
        return -1;
    }
}


int structurealign(int argc, const char **argv, const Command& command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    structureAlignDefault(par);
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);

    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);

    bool sameDB = false;
    uint16_t extended = DBReader<unsigned int>::getExtendedDbtype(FileUtil::parseDbType(par.db3.c_str()));
    bool alignmentIsExtended = extended & Parameters::DBTYPE_EXTENDED_INDEX_NEED_SRC;
    IndexReader tAADbr(par.db2, par.threads,
                             alignmentIsExtended ? IndexReader::SRC_SEQUENCES : IndexReader::SEQUENCES,
                             (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);

    std::string t3DiDbrName =  StructureUtil::getIndexWithSuffix(par.db2, "_ss");
    bool is3DiIdx = Parameters::isEqualDbtype(FileUtil::parseDbType(t3DiDbrName.c_str()),
                                              Parameters::DBTYPE_INDEX_DB);

    IndexReader t3DiDbr(is3DiIdx ? t3DiDbrName : par.db2, par.threads,
                              alignmentIsExtended ? IndexReader::SRC_SEQUENCES : IndexReader::SEQUENCES,
                              (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0,
                              DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA,
                              alignmentIsExtended ? "_seq_ss" : "_ss");

    IndexReader *q3DiDbr = NULL;
    IndexReader *qAADbr = NULL;

    if (par.db1.compare(par.db2) == 0) {
        sameDB = true;
        q3DiDbr = &t3DiDbr;
        qAADbr = &tAADbr;
    } else {
        qAADbr = new IndexReader(par.db1, par.threads, IndexReader::SRC_SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);
        q3DiDbr = new IndexReader(StructureUtil::getIndexWithSuffix(par.db1, "_ss"), par.threads, IndexReader::SRC_SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);
    }

    bool db1CaExist = FileUtil::fileExists((par.db1 + "_ca.dbtype").c_str());
    bool db2CaExist = FileUtil::fileExists((par.db2 + "_ca.dbtype").c_str());
    if(Parameters::isEqualDbtype(tAADbr.getDbtype(), Parameters::DBTYPE_INDEX_DB)){
        db2CaExist = true;
    }
    if(par.sortByStructureBits) {
        bool disableStructureBits = false;
        if(db1CaExist == false || db2CaExist == false){
            Debug(Debug::WARNING) << "Cannot find " << FileUtil::baseName(par.db1) << " C-alpha or " << FileUtil::baseName(par.db2) << " C-alpha database\n";
            disableStructureBits = true;
        }
        if(par.alignmentMode == 1 || par.alignmentMode == 2){
            Debug(Debug::WARNING) << "Cannot use --sort-by-structure-bits 1 with --alignment-mode 1 or 2\n";
            disableStructureBits = true;
        }
        if(disableStructureBits){
            Debug(Debug::WARNING) << "Disabling --sort-by-structure-bits\n";
            Debug(Debug::WARNING) << "This impacts the final score and ranking of hits, but not E-values themselves. Ranking alterations primarily occur for E-values < 10^-1.\n";
            par.sortByStructureBits = false;
        }
    }

    if (par.align2b && par.alignmentMode == Parameters::ALIGNMENT_MODE_SCORE_ONLY) {
        Debug(Debug::ERROR) << "--align2b needs the aligned region of the top hit, cannot use --alignment-mode 1\n";
        EXIT(EXIT_FAILURE);
    }
    if (par.align2bExclusion && par.align2b == false) {
        Debug(Debug::ERROR) << "--align2b-exclusion needs --align2b 1\n";
        EXIT(EXIT_FAILURE);
    }
    NcbiTaxonomy *taxonomy = NULL;
    MappingReader *qTaxMap = NULL;
    MappingReader *tTaxMap = NULL;
    if (par.align2bExclusion) {
        const std::string tDbName = PrefilteringIndexReader::dbPathWithoutIndex(par.db2);
        const std::string qDbName = PrefilteringIndexReader::dbPathWithoutIndex(par.db1);
        taxonomy = NcbiTaxonomy::openTaxonomy(tDbName);
        tTaxMap = new MappingReader(tDbName);
        qTaxMap = (qDbName == tDbName) ? tTaxMap : new MappingReader(qDbName);
    }

    DBReader<unsigned int> resultReader(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    resultReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    int dbtype =  Parameters::DBTYPE_ALIGNMENT_RES;
    if(alignmentIsExtended){
        dbtype = DBReader<unsigned int>::setExtendedDbtype(dbtype, Parameters::DBTYPE_EXTENDED_INDEX_NEED_SRC);
    }
    DBWriter dbw(par.db4.c_str(), par.db4Index.c_str(), static_cast<unsigned int>(par.threads), par.compressed,  dbtype);
    dbw.open();

    // align2b exclusion: per rank <db4>_<rank>, <db4>_<rank>_t1noexcl and <db4>_t1align_<rank>
    static const char *exclusionRankNames[] = {"species", "genus", "family", "order", "class", "phylum", "superkingdom"};
    const size_t exclusionRankCnt = par.align2bExclusion ? sizeof(exclusionRankNames) / sizeof(exclusionRankNames[0]) : 0;
    std::vector<int> exclusionRankIdx;
    std::vector<DBWriter*> rankDbw;
    std::vector<DBWriter*> rankT1NoExclDbw;
    std::vector<DBWriter*> rankT1AlignDbw;
    for (size_t r = 0; r < exclusionRankCnt; r++) {
        const std::string rank = exclusionRankNames[r];
        exclusionRankIdx.push_back(NcbiTaxonomy::findRankIndex(rank));
        const std::string names[3] = {par.db4 + "_" + rank, par.db4 + "_" + rank + "_t1noexcl", par.db4 + "_t1align_" + rank};
        std::vector<DBWriter*> *writers[3] = {&rankDbw, &rankT1NoExclDbw, &rankT1AlignDbw};
        for (size_t k = 0; k < 3; k++) {
            DBWriter *w = new DBWriter(names[k].c_str(), (names[k] + ".index").c_str(), static_cast<unsigned int>(par.threads), par.compressed, dbtype);
            w->open();
            writers[k]->push_back(w);
        }
    }

    bool needTMaligner = (par.tmScoreThr > 0);
    bool needLDDT = (par.lddtThr > 0);
    if (par.sortByStructureBits) {
        needLDDT = true;
        needTMaligner = true;
    } else {
        if (needTMaligner && (db1CaExist == false || db2CaExist == false)) {
            Debug(Debug::WARNING) << "Cannot use --tmscore-threshold with --sort-by-structure-bits 0\n"
                                  << "Disabling --tmscore-threshold\n";
            needTMaligner = false;
        }
        if (needLDDT && (db1CaExist == false || db2CaExist == false)) {
            Debug(Debug::WARNING) << "Cannot use --lddt-threshold with --sort-by-structure-bits 0\n"
                                  << "Disabling --lddt-threshold\n";
            needLDDT = false;
        }
    }
    bool needCalpha = (needTMaligner || needLDDT);
    IndexReader *qcadbr = NULL;
    IndexReader *tcadbr = NULL;
    if(needCalpha){
        qcadbr = new IndexReader(
                par.db1,
                par.threads,
                IndexReader::makeUserDatabaseType(LocalParameters::INDEX_DB_CA_KEY_DB1),
                touch ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0,
                DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA,
                "_ca");
        if (sameDB) {
            tcadbr = qcadbr;
        } else {
             tcadbr = new IndexReader(
                    par.db2,
                    par.threads,
                    alignmentIsExtended ? IndexReader::makeUserDatabaseType(LocalParameters::INDEX_DB_CA_KEY_DB2) :
                                           IndexReader::makeUserDatabaseType(LocalParameters::INDEX_DB_CA_KEY_DB1),
                    touch ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0,
                    DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA,
                    alignmentIsExtended ? "_seq_ca" : "_ca"
            );
        }
    }

    SubstitutionMatrix subMat3Di(par.scoringMatrixFile.values.aminoacid().c_str(), 2.1, par.scoreBias);
    std::string blosum;
    for (size_t i = 0; i < par.substitutionMatrices.size(); i++) {
        if (par.substitutionMatrices[i].name == "blosum62.out") {
            std::string matrixData((const char *)par.substitutionMatrices[i].subMatData, par.substitutionMatrices[i].subMatDataLen);
            std::string matrixName = par.substitutionMatrices[i].name;
            char * serializedMatrix = BaseMatrix::serialize(matrixName, matrixData);
            blosum.assign(serializedMatrix);
            free(serializedMatrix);
            break;
        }
    }
    float aaFactor = (par.alignmentType == LocalParameters::ALIGNMENT_TYPE_3DI_AA) ? 1.4 : 0.0;
    SubstitutionMatrix subMatAA(blosum.c_str(), aaFactor, par.scoreBias);
    //temporary output file
    Debug::Progress progress(resultReader.getSize());

    // sub. mat needed for query profile
    int8_t * tinySubMatAA = (int8_t*) mem_align(ALIGN_INT, subMatAA.alphabetSize * 32);
    int8_t * tinySubMat3Di = (int8_t*) mem_align(ALIGN_INT, subMat3Di.alphabetSize * 32);

    for (int i = 0; i < subMat3Di.alphabetSize; i++) {
        for (int j = 0; j < subMat3Di.alphabetSize; j++) {
            tinySubMat3Di[i * subMat3Di.alphabetSize + j] = subMat3Di.subMatrix[i][j]; // for farrar profile
        }
    }
    for (int i = 0; i < subMatAA.alphabetSize; i++) {
        for (int j = 0; j < subMatAA.alphabetSize; j++) {
            tinySubMatAA[i * subMatAA.alphabetSize + j] = subMatAA.subMatrix[i][j];
        }
    }

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        EvalueNeuralNet evaluer(tAADbr.sequenceReader->getAminoAcidDBSize(), &subMat3Di);
        std::vector<Matcher::result_t> alignmentResult;
        StructureSmithWaterman structureSmithWaterman(par.maxSeqLen, subMat3Di.alphabetSize, par.compBiasCorrection, par.compBiasCorrectionScale, &subMatAA, &subMat3Di);
        StructureSmithWaterman reverseStructureSmithWaterman(par.maxSeqLen, subMat3Di.alphabetSize, par.compBiasCorrection, par.compBiasCorrectionScale, &subMatAA, &subMat3Di);
        TMaligner *tmaligner = NULL;
        if(needTMaligner) {
            tmaligner = new TMaligner(
                    std::max(q3DiDbr->sequenceReader->getMaxSeqLen() + 1, t3DiDbr.sequenceReader->getMaxSeqLen() + 1), false, true, par.exactTMscore);
        }
        LDDTCalculator *lddtcalculator = NULL;
        if(needLDDT) {
            lddtcalculator = new LDDTCalculator(q3DiDbr->sequenceReader->getMaxSeqLen() + 1,  t3DiDbr.sequenceReader->getMaxSeqLen() + 1);
        }
        Sequence qSeqAA(par.maxSeqLen, qAADbr->getDbtype(), (const BaseMatrix *) &subMatAA, 0, false, par.compBiasCorrection);
        Sequence qSeq3Di(par.maxSeqLen, q3DiDbr->getDbtype(), (const BaseMatrix *) &subMat3Di, 0, false, par.compBiasCorrection);
        Sequence tSeqAA(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, (const BaseMatrix *) &subMatAA, 0, false, par.compBiasCorrection);
        Sequence tSeq3Di(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, (const BaseMatrix *) &subMat3Di, 0, false, par.compBiasCorrection);
        Sequence t1SeqAA(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, (const BaseMatrix *) &subMatAA, 0, false, par.compBiasCorrection);
        Sequence t1Seq3Di(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, (const BaseMatrix *) &subMat3Di, 0, false, par.compBiasCorrection);
        // align2b: primary pass-1 hits (score = raw bitscore), the candidate list with the rank of
        // each candidate's LCA with the query, and one T1-region rescan per distinct T1
        struct Rescan {
            size_t t1Idx;
            size_t scanned;
            std::vector<std::pair<size_t, Matcher::result_t>> passed;
        };
        std::vector<Matcher::result_t> p1Hits;
        std::vector<unsigned int> candKeys;
        std::vector<int> candLcaRank;
        std::vector<Rescan> rescans;
        std::vector<Matcher::result_t> rankResults;
        std::vector<std::string> rankBuffers(3 * exclusionRankCnt);
        auto sortHits = [&](std::vector<Matcher::result_t> &hits) {
            if (hits.size() > 1) {
                if(par.sortByStructureBits) {
                    SORT_SERIAL(hits.begin(), hits.end(), compareHitsByStructureBits);
                } else {
                    SORT_SERIAL(hits.begin(), hits.end(), Matcher::compareHits);
                }
            }
        };
        std::string backtrace;
        char buffer[1024+32768];
        std::string resultBuffer;

        Coordinate16 qcoords;
        Coordinate16 tcoords;

        TMaligner::TMscoreResult tmres;
        LDDTCalculator::LDDTScoreResult lddtres;
        // write output file

#pragma omp for schedule(dynamic, 1)
        for (size_t id = 0; id < resultReader.getSize(); id++) {
            progress.updateProgress();
            char *data = resultReader.getData(id, thread_idx);
            char *origData = data;
            size_t queryKey = resultReader.getDbKey(id);
            if(*data != '\0') {
                unsigned int queryId = q3DiDbr->sequenceReader->getId(queryKey);

                char *querySeqAA = qAADbr->sequenceReader->getData(queryId, thread_idx);
                char *querySeq3Di = q3DiDbr->sequenceReader->getData(queryId, thread_idx);
                unsigned int querySeqLen = q3DiDbr->sequenceReader->getSeqLen(queryId);
                qSeq3Di.mapSequence(id, queryKey, querySeq3Di, querySeqLen);
                qSeqAA.mapSequence(id, queryKey, querySeqAA, querySeqLen);
                if(needCalpha){
                    size_t qId = qcadbr->sequenceReader->getId(queryKey);
                    char *qcadata = qcadbr->sequenceReader->getData(qId, thread_idx);
                    size_t qCaLength = qcadbr->sequenceReader->getEntryLen(qId);
                    float* queryCaData = qcoords.read(qcadata, qSeq3Di.L, qCaLength);
                    if(needTMaligner){
                        tmaligner->initQuery(queryCaData, &queryCaData[qSeq3Di.L], &queryCaData[qSeq3Di.L+qSeq3Di.L], NULL, qSeq3Di.L);
                    }
                    if(needLDDT){
                        lddtcalculator->initQuery(qSeq3Di.L, queryCaData, &queryCaData[qSeq3Di.L], &queryCaData[qSeq3Di.L+qSeq3Di.L]);
                    }
                }
                std::pair<double, double> muLambda = evaluer.predictMuLambda(qSeq3Di.numSequence, qSeq3Di.L);
                structureSmithWaterman.ssw_init(&qSeqAA, &qSeq3Di, tinySubMatAA, tinySubMat3Di, &subMatAA);
                qSeq3Di.reverse();
                qSeqAA.reverse();
                reverseStructureSmithWaterman.ssw_init(&qSeqAA, &qSeq3Di, tinySubMatAA, tinySubMat3Di, &subMatAA);
                int passedNum = 0;
                int rejected = 0;
                while (*data != '\0' && passedNum < par.maxAccept && rejected < par.maxRejected) {
                    char dbKeyBuffer[255 + 1];
                    Util::parseKey(data, dbKeyBuffer);
                    data = Util::skipLine(data);
                    const unsigned int dbKey = (unsigned int) strtoul(dbKeyBuffer, NULL, 10);
                    unsigned int targetId = t3DiDbr.sequenceReader->getId(dbKey);
                    const bool isIdentity = (queryId == targetId && (par.includeIdentity || sameDB))? true : false;

                    char * targetSeq3Di = t3DiDbr.sequenceReader->getData(targetId, thread_idx);
                    char * targetSeqAA = tAADbr.sequenceReader->getData(targetId, thread_idx);
                    const int targetSeqLen = static_cast<int>(t3DiDbr.sequenceReader->getSeqLen(targetId));

                    tSeq3Di.mapSequence(targetId, dbKey, targetSeq3Di, targetSeqLen);
                    tSeqAA.mapSequence(targetId, dbKey, targetSeqAA, targetSeqLen);
                    if(Util::canBeCovered(par.covThr, par.covMode, qSeq3Di.L, targetSeqLen) == false){
                        rejected++;
                        continue;
                    }
                    Matcher::result_t res;
                    if(alignStructure(structureSmithWaterman, reverseStructureSmithWaterman,
                                      tSeqAA, tSeq3Di, querySeqLen, targetSeqLen,
                                      evaluer, muLambda, res, backtrace, par,
                                      par.evalThr, par.covThr, INT32_MIN) == -1){
                        rejected++;
                        continue;
                    }

                    if (Alignment::checkCriteria(res, isIdentity, par.evalThr, par.seqIdThr, par.alnLenThr, par.covMode, par.covThr)) {
                        // raw bitscore, before it gets scaled into structure bits
                        const int32_t rawScore = res.score;
                        if(needCalpha) {
                            size_t tId = tcadbr->sequenceReader->getId(res.dbKey);
                            char *tcadata = tcadbr->sequenceReader->getData(tId, thread_idx);
                            size_t tCaLength = tcadbr->sequenceReader->getEntryLen(tId);
                            float* targetCaData = tcoords.read(tcadata, res.dbLen, tCaLength);
                            if(needTMaligner) {
                                tmres = tmaligner->computeTMscore(targetCaData,
                                                                  &targetCaData[res.dbLen],
                                                                  &targetCaData[res.dbLen +
                                                                                res.dbLen],
                                                                  res.dbLen,
                                                                  res.qStartPos,
                                                                  res.dbStartPos,
                                                                  res.backtrace,
                                                                  TMaligner::normalization(par.tmScoreThrMode, std::min(res.qEndPos - res.qStartPos, res.dbEndPos - res.dbStartPos ), res.qLen, res.dbLen));

                                if (tmres.tmscore < par.tmScoreThr) {
                                    continue;
                                }
                            }
                            if(needLDDT){
                                lddtres = lddtcalculator->computeLDDTScore(res.dbLen, res.qStartPos, res.dbStartPos,
                                                                           res.backtrace,
                                                                           targetCaData, &targetCaData[res.dbLen],
                                                                           &targetCaData[res.dbLen+res.dbLen]);

                                if(lddtres.avgLddtScore < par.lddtThr){
                                    continue;
                                }
                                res.dbcov = lddtres.avgLddtScore;
                            }
                            if(par.sortByStructureBits && needTMaligner && needLDDT){
                                res.score = res.score * sqrt(lddtres.avgLddtScore * tmres.tmscore);
                            }
                        }


                        if (par.align2b) {
                            p1Hits.emplace_back(res);
                            p1Hits.back().score = rawScore;
                        }
                        alignmentResult.emplace_back(res);
                        int altAli = par.altAlignment;
                        bool moreAltAli = true;
                        while(altAli && moreAltAli){
                            Matcher::result_t altRes;
                            if(computeAlternativeAlignment(structureSmithWaterman, reverseStructureSmithWaterman,
                                                           tSeqAA, tSeq3Di, querySeqLen, targetSeqLen,
                                                           evaluer, muLambda, res, altRes,
                                                           backtrace, par) == -1) {
                                moreAltAli = false;
                                continue;
                            }
                            alignmentResult.push_back(altRes);
                            res = altRes;
                            altAli--;
                        }
                        passedNum++;
                        rejected = 0;
                    } else {
                        rejected++;
                    }
                }

                // 2b: realign the aligned region of the top hit (T1) against all candidates and
                // keep targets with bits(T1,T) >= bits(Q,T1). Reported alignments are T1 vs target.
                if (par.align2b && p1Hits.empty() == false) {
                    // exclusion can't work without taxonomy, so a missing taxon is an error
                    auto lookupTaxon = [&](MappingReader *mapping, size_t key, const char *what) -> TaxID {
                        const TaxID taxon = mapping->lookup(key);
                        if (taxon == 0 || taxonomy->nodeExists(taxon) == false) {
                            Debug(Debug::ERROR) << what << " " << key << " has no taxonomy assigned, cannot use --align2b-exclusion\n";
                            EXIT(EXIT_FAILURE);
                        }
                        return taxon;
                    };
                    TaxID queryTax = 0;
                    if (par.align2bExclusion) {
                        queryTax = lookupTaxon(qTaxMap, queryKey, "Query");
                    }

                    // candidates = the query's prefilter list; their LCA rank with the query doesn't
                    // depend on T1 or the exclusion rank, so it is computed once here
                    data = origData;
                    while (*data != '\0') {
                        char dbKeyBuffer[255 + 1];
                        Util::parseKey(data, dbKeyBuffer);
                        data = Util::skipLine(data);
                        const unsigned int dbKey = (unsigned int) strtoul(dbKeyBuffer, NULL, 10);
                        candKeys.push_back(dbKey);
                        candLcaRank.push_back(par.align2bExclusion
                                              ? lcaRankIndex(taxonomy, queryTax, lookupTaxon(tTaxMap, dbKey, "Target"))
                                              : -1);
                    }

                    // T1 = best raw bitscore among the pass-1 hits, optionally skipping targets within rankIdx of the query
                    auto pickT1 = [&](int rankIdx) -> int {
                        int best = -1;
                        for (size_t i = 0; i < p1Hits.size(); i++) {
                            if (rankIdx >= 0 && isTooClose(lcaRankIndex(taxonomy, queryTax, lookupTaxon(tTaxMap, p1Hits[i].dbKey, "Target")), rankIdx)) {
                                continue;
                            }
                            if (best < 0 || p1Hits[i].score > p1Hits[best].score) {
                                best = static_cast<int>(i);
                            }
                        }
                        return best;
                    };

                    // with stopEarly the max-accept/max-rejected limits apply while scanning, otherwise
                    // all candidates are scanned once and select() applies the limits per output
                    auto rescan = [&](size_t t1Idx, bool stopEarly) -> size_t {
                        for (size_t i = 0; i < rescans.size(); i++) {
                            if (rescans[i].t1Idx == t1Idx) {
                                return i;
                            }
                        }
                        rescans.emplace_back();
                        Rescan &scan = rescans.back();
                        scan.t1Idx = t1Idx;
                        scan.scanned = 0;

                        const Matcher::result_t &t1Hit = p1Hits[t1Idx];
                        unsigned int t1Id = t3DiDbr.sequenceReader->getId(t1Hit.dbKey);
                        char *t1DataAA = tAADbr.sequenceReader->getData(t1Id, thread_idx);
                        char *t1Data3Di = t3DiDbr.sequenceReader->getData(t1Id, thread_idx);
                        const unsigned int t1RegionLen = t1Hit.dbEndPos - t1Hit.dbStartPos + 1;
                        t1SeqAA.mapSequence(t1Id, t1Hit.dbKey, t1DataAA + t1Hit.dbStartPos, t1RegionLen);
                        t1Seq3Di.mapSequence(t1Id, t1Hit.dbKey, t1Data3Di + t1Hit.dbStartPos, t1RegionLen);
                        std::pair<double, double> t1MuLambda = evaluer.predictMuLambda(t1Seq3Di.numSequence, t1Seq3Di.L);
                        structureSmithWaterman.ssw_init(&t1SeqAA, &t1Seq3Di, tinySubMatAA, tinySubMat3Di, &subMatAA);
                        t1Seq3Di.reverse();
                        t1SeqAA.reverse();
                        reverseStructureSmithWaterman.ssw_init(&t1SeqAA, &t1Seq3Di, tinySubMatAA, tinySubMat3Di, &subMatAA);

                        int scanPassed = 0;
                        int scanRejected = 0;
                        for (size_t c = 0; c < candKeys.size(); c++) {
                            if (stopEarly && (scanPassed >= par.maxAccept || scanRejected >= par.maxRejected)) {
                                break;
                            }
                            scan.scanned = c + 1;
                            const unsigned int dbKey = candKeys[c];
                            unsigned int targetId = t3DiDbr.sequenceReader->getId(dbKey);
                            char * targetSeq3Di = t3DiDbr.sequenceReader->getData(targetId, thread_idx);
                            char * targetSeqAA = tAADbr.sequenceReader->getData(targetId, thread_idx);
                            const int targetSeqLen = static_cast<int>(t3DiDbr.sequenceReader->getSeqLen(targetId));
                            tSeq3Di.mapSequence(targetId, dbKey, targetSeq3Di, targetSeqLen);
                            tSeqAA.mapSequence(targetId, dbKey, targetSeqAA, targetSeqLen);

                            // cut on bitscore only
                            Matcher::result_t res;
                            if (alignStructure(structureSmithWaterman, reverseStructureSmithWaterman,
                                               tSeqAA, tSeq3Di, t1RegionLen, targetSeqLen,
                                               evaluer, t1MuLambda, res, backtrace, par,
                                               DBL_MAX, 0.0f, t1Hit.score) == -1) {
                                scanRejected++;
                                continue;
                            }
                            scan.passed.emplace_back(c, res);
                            scanPassed++;
                            scanRejected = 0;
                        }
                        return rescans.size() - 1;
                    };

                    // exclusion (rankIdx < 0: none) plus the max-accept/max-rejected limits, as if
                    // excluded candidates were not in the candidate list
                    auto select = [&](size_t scanIdx, int rankIdx, std::vector<Matcher::result_t> &out) {
                        const Rescan &scan = rescans[scanIdx];
                        out.clear();
                        int selPassed = 0;
                        int selRejected = 0;
                        size_t p = 0;
                        for (size_t c = 0; c < scan.scanned && selPassed < par.maxAccept && selRejected < par.maxRejected; c++) {
                            const bool passed = (p < scan.passed.size() && scan.passed[p].first == c);
                            const bool excluded = (rankIdx >= 0 && isTooClose(candLcaRank[c], rankIdx));
                            if (excluded == false) {
                                if (passed) {
                                    out.push_back(scan.passed[p].second);
                                    selPassed++;
                                    selRejected = 0;
                                } else {
                                    selRejected++;
                                }
                            }
                            if (passed) {
                                p++;
                            }
                        }
                    };

                    // plain align2b; t1noexcl outputs reuse its rescan
                    const int baseT1 = pickT1(-1);
                    const size_t baseScan = rescan(baseT1, par.align2bExclusion == false);
                    select(baseScan, -1, alignmentResult);

                    if (par.align2bExclusion) {
                        for (size_t r = 0; r < exclusionRankCnt; r++) {
                            const int rankIdx = exclusionRankIdx[r];
                            const int exclT1 = pickT1(rankIdx);
                            if (exclT1 >= 0) {
                                select(rescan(exclT1, false), rankIdx, rankResults);
                                sortHits(rankResults);
                                for (size_t i = 0; i < rankResults.size(); i++) {
                                    size_t len = Matcher::resultToBuffer(buffer, rankResults[i], par.addBacktrace);
                                    rankBuffers[3 * r].append(buffer, len);
                                }
                                size_t len = Matcher::resultToBuffer(buffer, p1Hits[exclT1], par.addBacktrace);
                                rankBuffers[3 * r + 2].append(buffer, len);
                            }
                            select(baseScan, rankIdx, rankResults);
                            sortHits(rankResults);
                            for (size_t i = 0; i < rankResults.size(); i++) {
                                size_t len = Matcher::resultToBuffer(buffer, rankResults[i], par.addBacktrace);
                                rankBuffers[3 * r + 1].append(buffer, len);
                            }
                        }
                    }
                }
            }


            sortHits(alignmentResult);
            for (size_t result = 0; result < alignmentResult.size(); result++) {
                size_t len = Matcher::resultToBuffer(buffer, alignmentResult[result], par.addBacktrace);
                resultBuffer.append(buffer, len);
            }
            dbw.writeData(resultBuffer.c_str(), resultBuffer.length(), queryKey, thread_idx);
            resultBuffer.clear();
            alignmentResult.clear();
            for (size_t r = 0; r < exclusionRankCnt; r++) {
                rankDbw[r]->writeData(rankBuffers[3 * r].c_str(), rankBuffers[3 * r].length(), queryKey, thread_idx);
                rankT1NoExclDbw[r]->writeData(rankBuffers[3 * r + 1].c_str(), rankBuffers[3 * r + 1].length(), queryKey, thread_idx);
                rankT1AlignDbw[r]->writeData(rankBuffers[3 * r + 2].c_str(), rankBuffers[3 * r + 2].length(), queryKey, thread_idx);
            }
            for (size_t i = 0; i < rankBuffers.size(); i++) {
                rankBuffers[i].clear();
            }
            p1Hits.clear();
            candKeys.clear();
            candLcaRank.clear();
            rescans.clear();
        }
        if(needTMaligner){
            delete tmaligner;
        }
        if(needLDDT){
            delete lddtcalculator;
        }
    }

    free(tinySubMatAA);
    free(tinySubMat3Di);

    dbw.close();
    resultReader.close();

    for (size_t r = 0; r < exclusionRankCnt; r++) {
        rankDbw[r]->close();
        rankT1NoExclDbw[r]->close();
        rankT1AlignDbw[r]->close();
        delete rankDbw[r];
        delete rankT1NoExclDbw[r];
        delete rankT1AlignDbw[r];
    }
    if (par.align2bExclusion) {
        if (qTaxMap != tTaxMap) {
            delete qTaxMap;
        }
        delete tTaxMap;
        delete taxonomy;
    }

    if(needCalpha){
        if (sameDB == false) {
            delete tcadbr;
        }
        delete qcadbr;
    }

    if (sameDB == false) {
        delete q3DiDbr;
        delete qAADbr;
    }

    return EXIT_SUCCESS;
}
