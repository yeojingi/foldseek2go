#ifndef FUNC_READER_H
#define FUNC_READER_H

#include "Util.h"
#include "MemoryMapped.h"
#include <algorithm>

#include <iostream>

class FuncReader {
public:
    // static std::pair<char *, size_t> serialize(const FuncReader &reader) {
    //     size_t serialized_size = reader.magicLen + reader.count * sizeof(Pair);
    //     char* data = (char*)malloc(serialized_size);
    //     memcpy(data, reader.magic, reader.magicLen);
    //     memcpy(data + reader.magicLen, reader.entries, reader.count * sizeof(Pair));
    //     return std::make_pair(data, serialized_size);
    // }

    FuncReader(const std::string &db) {
        std::string input = db + "_func";
        file = new MemoryMapped(input, MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
        if (!file->isValid()) {
            delete file;
            file = NULL;
            Debug(Debug::ERROR) << db << "_func does not exist. Please create the function mapping with createfuncdb!\n";
            EXIT(EXIT_FAILURE);
        }

        char *data = (char *) file->getData();
        size_t dataSize = file->size();
        // if (file->size() > magicLen && memcmp(data, magic, magicLen) == 0) {
        //     entries = reinterpret_cast<Pair*>(data + magicLen);
        //     count = (dataSize - magicLen) / sizeof(Pair);
        //     return;
        // }
        
        std::vector<std::pair<unsigned int, unsigned int>> mapping;
        size_t currPos = 0;
        const char *cols[3];
        size_t isSorted = true;
        unsigned int prevId = 0;
        while (currPos < dataSize) {
            Util::getWordsOfLine(data, cols, 2);
            unsigned int id = Util::fast_atoi<size_t>(cols[0]);
            isSorted *= (id >= prevId);
            unsigned int goid = Util::fast_atoi<size_t>(cols[1]);
            data = Util::skipLine(data);
            mapping.push_back(std::make_pair(id, goid));
            currPos = data - (char *) file->getData();
            prevId = id;
            // std::cout << id << " " << goid << std::endl;
        }
        file->close();
        delete file;
        file = NULL;
        if (mapping.size() == 0) {
            Debug(Debug::ERROR) << db << "_func is empty. Rerun createtaxdb to recreate function mapping with createfuncdb.\n";
            EXIT(EXIT_FAILURE);
        }

        count = mapping.size();
        entries = new Pair[count];
        unsigned int prev_id = 0;
        entries[0].dbkey = mapping[0].first;
        for (size_t i = 0; i < count; i++) {
            if (entries[prev_id].dbkey != mapping[i].first) {
                entries[i].dbkey = mapping[i].first;
                entries[i].goids.push_back(mapping[i].second);
                prev_id = i;
            } else {
                entries[prev_id].goids.push_back(mapping[i].second);
                // std::cout << mapping[i].second << " " << entries[prev_id].goids.size() << std::endl;
            }
        }
        if (isSorted == false) {
            std::stable_sort(entries, entries + count, compareTaxa);
        }

    }

    ~FuncReader() {
        if (file != NULL) {
            file->close();
            delete file;
        } else {
            delete[] entries;
        }
    }

    std::vector<unsigned int> * lookup(unsigned int key) {
        std::vector<unsigned int> *goids;
        // match dbKey to its taxon based on mapping
        Pair val;
        val.dbkey = key;
        Pair* end = entries + count;
        Pair* found = std::upper_bound(entries, end, val, compareTaxa);
        if (found == end || found->dbkey != key) {
            goids = 0;
        } else {
            goids = &(found->goids);
            // std::cout << key << "-" << *goids[0] << std::endl;
        }

        // std::cout << entries[key].dbkey << " " << entries[key].goids[0] << std::endl;

        // goids = &entries[key].goids;

        return goids;
    }

private:
    MemoryMapped* file;
    struct __attribute__((__packed__)) Pair{
        unsigned int dbkey;
        // unsigned int taxon;
        std::vector<unsigned int> goids;
    };
    Pair* entries;
    size_t count;
    //                    T  A   X   M  Version
    // const char magic[5] = {19, 0, 23, 12, 0};
    // const size_t magicLen = 5;
    static bool compareTaxa(const Pair &lhs, const Pair &rhs) {
        return (lhs.dbkey <= rhs.dbkey);
    }
};

#endif
