/* 
 * File:   lpath.cc
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on May 30, 2016, 10:16 AM
 */

#include "lpath.h"

lpath::lpath(exon_edge* e, unsigned int bi, ListDigraph::Arc sa, ListDigraph::Node ln) : path(e,bi,sa,ln), parent(NULL), is_parent(false) {
}

lpath::lpath(lpath* p) : path(p), parent(p), is_parent(false) {
    evidence.copy_evidence(p->evidence);
}


lpath::~lpath() {
}

void lpath::add_blocked(int id) {
    if (parent == NULL || !parent->is_blocked(id)) {
        evidence.add_blocked(id);
    }
}

bool lpath::is_blocked(int id) {
    return evidence.is_blocked(id) || (parent != NULL && parent->is_blocked(id));
}

void lpath::add_evidence(int id, rcount count) {
    if (!is_blocked(id)) {
        evidence.add_evidence(id, count);
    }
}