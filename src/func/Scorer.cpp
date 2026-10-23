#include "Scorer.h"
#include "Debug.h"
#include "Matcher.h"
#include "IndexReader.h"
#include "BaseMatrix.h"
#include <algorithm>

void ScoreAlignments(DBReader<unsigned int>* alnDbr, EvidenceScore* evidenceScore, IndexReader* tGoDbr) {
    const int thread_idx = 0;
    Debug::Progress progress(alnDbr->getSize());

    for (size_t i = 0; i < alnDbr->getSize(); i++) {
        progress.updateProgress();
        const unsigned int queryKey = alnDbr->getDbKey(i);

        char *data = alnDbr->getData(i, thread_idx);

        while (*data != '\0') {
            Matcher::result_t res = Matcher::parseAlignmentRecord(data, true);
            data = Util::skipLine(data);
            size_t tGoId = tGoDbr->sequenceReader->getId(res.dbKey);
            if (tGoId == UINT_MAX) continue; // GO annotation 없는 target 제외
            evidenceScore->addScore(queryKey, tGoId, res.eval);
        }
    }
}

void ScoreAlignmentsBitscore(DBReader<unsigned int>* alnDbr, EvidenceScore* bsScore, IndexReader* tGoDbr) {
    const int thread_idx = 0;
    Debug::Progress progress(alnDbr->getSize());
    for (size_t i = 0; i < alnDbr->getSize(); i++) {
        progress.updateProgress();
        const unsigned int queryKey = alnDbr->getDbKey(i);
        char *data = alnDbr->getData(i, thread_idx);
        while (*data != '\0') {
            Matcher::result_t res = Matcher::parseAlignmentRecord(data, true);
            data = Util::skipLine(data);
            size_t tGoId = tGoDbr->sequenceReader->getId(res.dbKey);
            if (tGoId == UINT_MAX) continue;
            bsScore->addScore(queryKey, tGoId, (float)res.score);
        }
    }
}

void ScoreAlignmentsSimilarity(DBReader<unsigned int>* alnDbr, EvidenceScore* simScore, IndexReader* tGoDbr, IndexReader* qDbr, IndexReader* tDbr, BaseMatrix* subMat) {
    const int thread_idx = 0;
    Debug::Progress progress(alnDbr->getSize());

    for (size_t i = 0; i < alnDbr->getSize(); i++) {
        progress.updateProgress();
        const unsigned int queryKey = alnDbr->getDbKey(i);

        char *data = alnDbr->getData(i, thread_idx);

        while (*data != '\0') {
            Matcher::result_t res = Matcher::parseAlignmentRecord(data, false);
            data = Util::skipLine(data);

            if (res.backtrace.empty()) {
                Debug(Debug::ERROR) << "Alignment result is missing backtrace information. Please rerun mmseqs search/align with the -a option.\n";
                EXIT(EXIT_FAILURE);
            }

            size_t tGoId = tGoDbr->sequenceReader->getId(res.dbKey);
            if (tGoId == UINT_MAX) continue; // GO annotation 없는 target 제외

            size_t qIdx = qDbr->sequenceReader->getId(queryKey);
            size_t tIdx = tDbr->sequenceReader->getId(res.dbKey);
            const char* qSeq = qDbr->sequenceReader->getData(qIdx, thread_idx);
            const char* tSeq = tDbr->sequenceReader->getData(tIdx, thread_idx);

            size_t qPos = (size_t)std::max(res.qStartPos, 0);
            size_t tPos = (size_t)std::max(res.dbStartPos, 0);
            int positives = 0;
            int total = 0;
            for (char op : res.backtrace) {
                total++;
                if (op == 'M') {
                    unsigned char qAA = subMat->aa2num[(unsigned char)qSeq[qPos]];
                    unsigned char tAA = subMat->aa2num[(unsigned char)tSeq[tPos]];
                    if (subMat->subMatrix[qAA][tAA] > 0) positives++;
                    qPos++;
                    tPos++;
                } else if (op == 'I') {
                    qPos++;
                } else if (op == 'D') {
                    tPos++;
                }
            }

            float similarity = (total > 0) ? (100.0f * (float)positives / (float)total) : 0.0f;
            simScore->addScore(queryKey, tGoId, similarity);
        }
    }
}