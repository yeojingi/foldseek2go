#include <climits>
#include <cmath>
#include <cstdio>

#include "LocalParameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "Coordinate16.h"

#ifdef OPENMP
#include <omp.h>
#endif

// Globularity from C-alpha coordinates, adapted from cath-alphaflow measure-globularity:
// packing density = mean number of residues with a C-alpha within --glob-distance of a
// hydrophobic residue's C-alpha (itself included), normed radius of gyration = C-alpha Rg
// divided by the Rg of a solid sphere with the volume of the protein.

// average volume of a residue in a folded protein (Angstrom^3), used for the protein volume
static const double RESIDUE_VOLUME = 135.0;

// same set as cath-alphaflow
static bool isHydrophobic(char aa) {
    switch (aa & ~0x20) {
        case 'G': case 'A': case 'L': case 'I': case 'V': case 'P': case 'F': case 'M': case 'W':
            return true;
        default:
            return false;
    }
}

int computeglob(int argc, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> db(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    db.open(DBReader<unsigned int>::NOSORT);

    std::string dbCa = par.db1 + "_ca";
    std::string dbCaIndex = par.db1 + "_ca.index";
    DBReader<unsigned int> dbCaReader(dbCa.c_str(), dbCaIndex.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    dbCaReader.open(DBReader<unsigned int>::NOSORT);

    DBWriter dbw(par.db2.c_str(), par.db2Index.c_str(), static_cast<unsigned int>(par.threads), par.compressed, Parameters::DBTYPE_GENERIC_DB);
    dbw.open();

    const float distSq = par.globDistance * par.globDistance;
    Debug::Progress progress(db.getSize());
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        Coordinate16 coords;
        char buffer[256];

#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < db.getSize(); ++i) {
            progress.updateProgress();
            const unsigned int key = db.getDbKey(i);
            const char* seqData = db.getData(i, thread_idx);
            const size_t seqLen = std::max(db.getEntryLen(i), (size_t)2) - 2;

            unsigned int caId = dbCaReader.getId(key);
            if (caId == UINT_MAX) {
                Debug(Debug::ERROR) << "Entry " << key << " has no C-alpha coordinates in " << dbCa << "\n";
                EXIT(EXIT_FAILURE);
            }
            const char* caData = dbCaReader.getData(caId, thread_idx);
            const size_t caLen = dbCaReader.getEntryLen(caId);
            float* ca = coords.read(caData, seqLen, caLen);
            const float* x = ca;
            const float* y = ca + seqLen;
            const float* z = ca + 2 * seqLen;

            // packing density
            size_t hydrophobicCnt = 0;
            size_t neighborSum = 0;
            for (size_t r = 0; r < seqLen; ++r) {
                if (isHydrophobic(seqData[r]) == false) {
                    continue;
                }
                hydrophobicCnt++;
                for (size_t s = 0; s < seqLen; ++s) {
                    const float dx = x[r] - x[s];
                    const float dy = y[r] - y[s];
                    const float dz = z[r] - z[s];
                    if (dx * dx + dy * dy + dz * dz <= distSq) {
                        neighborSum++;
                    }
                }
            }
            const double packingDensity = (hydrophobicCnt > 0)
                                          ? static_cast<double>(neighborSum) / hydrophobicCnt : NAN;

            // radius of gyration, all C-alphas weighted equally
            double cx = 0.0, cy = 0.0, cz = 0.0;
            for (size_t r = 0; r < seqLen; ++r) {
                cx += x[r];
                cy += y[r];
                cz += z[r];
            }
            double normedRg = NAN;
            if (seqLen > 0) {
                cx /= seqLen;
                cy /= seqLen;
                cz /= seqLen;
                double sumSq = 0.0;
                for (size_t r = 0; r < seqLen; ++r) {
                    const double dx = x[r] - cx;
                    const double dy = y[r] - cy;
                    const double dz = z[r] - cz;
                    sumSq += dx * dx + dy * dy + dz * dz;
                }
                const double rg = std::sqrt(sumSq / seqLen);
                // solid sphere with the same volume: r = (3V / 4pi)^(1/3), Rg = sqrt(3/5) r
                const double radius = std::cbrt(3.0 * RESIDUE_VOLUME * seqLen / (4.0 * M_PI));
                normedRg = rg / (std::sqrt(3.0 / 5.0) * radius);
            }

            const int len = snprintf(buffer, sizeof(buffer), "%.3f\t%.3f\n", packingDensity, normedRg);
            dbw.writeData(buffer, len, key, thread_idx);
        }
    }

    dbw.close();
    dbCaReader.close();
    db.close();

    return EXIT_SUCCESS;
}
