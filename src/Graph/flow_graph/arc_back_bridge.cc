/* 
 * File:   arc_back_bridge.cc
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on June 15, 2016, 2:59 PM
 */

#include "arc_back_bridge.h"

arc_back_bridge::arc_back_bridge() {
}

arc_back_bridge::~arc_back_bridge() {
}

bool arc_back_bridge::is_evidenced_path() {
    return !bridges->empty();
}

bool arc_back_bridge::has_path_evidence(int id) {
    return bridges->find(id) != bridges->end();
}


arc_back_bridge::iterator arc_back_bridge::begin() {
    return bridges->begin();
}

arc_back_bridge::iterator arc_back_bridge::end() {
    return bridges->end();
}

int arc_back_bridge::size() {
    return bridges->size();
}

void arc_back_bridge::clear() {
    return bridges->clear();
}

void arc_back_bridge::remove_id(int id) {
    bridges->erase(id);
}

void arc_back_bridge::replace_evidence(int old_id, int new_id) {
    path_evidence_set<int>::iterator it = bridges->find(old_id);
    if (it == bridges->end()) {
        return; // nothing to do here
    }
    
    bridges->insert(new_id);
    bridges->erase(it);
}

void arc_back_bridge::add_evidence_if(int old_id, int new_id) {
 path_evidence_set<int>::iterator it = bridges->find(old_id);
    if (it == bridges->end()) {
        return; // nothing to do here
    }
    
    bridges->insert(new_id);
}