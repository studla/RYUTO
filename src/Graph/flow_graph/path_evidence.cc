/* 
 * File:   path_evidence.cc
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on May 24, 2016, 10:41 AM
 */

#include "path_evidence.h"

path_evidence::path_evidence()  {
}

path_evidence::~path_evidence() {
}

void path_evidence::add_evidence(int id, rcount count) {
    // this should never be called without checking all blocks first!
    evidence[id] += count;
}

void path_evidence::add_blocked(int id) {
    blocked.insert(id);
    evidence.erase(id);
}

bool path_evidence::is_blocked(int id) {
    return blocked.find(id) != blocked.end();
}

bool path_evidence::is_evidence(int id) {
    return evidence.find(id) != evidence.end();
}

void path_evidence::copy_evidence(path_evidence& other) {
    evidence = other.evidence; // deep copies
}

path_evidence_map<int, rcount>::iterator path_evidence::begin() {
    return evidence.begin();
}

path_evidence_map<int, rcount>::iterator path_evidence::end() {
    return evidence.end();
}

bool path_evidence::no_evidence() {
    return evidence.empty();
}