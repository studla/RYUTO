/* 
 * File:   arc_bridge.cc
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on June 7, 2016, 8:56 AM
 */

#include "arc_bridge.h"
#include "arc_back_bridge.h"

arc_bridge::arc_bridge() {
}

arc_bridge::~arc_bridge() {
}

bool arc_bridge::is_evidenced_path() {
    return !bridges->empty();
}

bool arc_bridge::has_path_evidence(int id) {
    return bridges->find(id) != bridges->end();
}

rcount arc_bridge::get_path_evidence(int id) {
   
    path_evidence_map<int, rcount>::iterator it = bridges->find(id);
    if (it != bridges->end()) {
        return it->second;
    }
    
    return 0;
}

void arc_bridge::replace_evidence(int old_id, int new_id) {
    path_evidence_map<int, rcount>::iterator it = bridges->find(old_id);
    if (it == bridges->end()) {
        return; // nothing to do here
    }
    
    bridges.ref()[new_id] += it->second;
    bridges->erase(it);
}

void arc_bridge::add_evidence_if(int old_id, int new_id) {
    path_evidence_map<int, rcount>::iterator it = bridges->find(old_id);
    if (it == bridges->end()) {
        return; // nothing to do here
    }
    
    bridges.ref()[new_id] += it->second; 
    it->second = 1;     // arbitrary evidence value for left over arcs
}

arc_bridge::iterator arc_bridge::begin() {
    return bridges->begin();
}

arc_bridge::iterator arc_bridge::end() {
    return bridges->end();
}


int arc_bridge::size() {
    return bridges->size();
}

void arc_bridge::clear() {
    return bridges->clear();
}

void arc_bridge::remove_id(int id) {
    bridges->erase(id);
}

rcount& arc_bridge::operator [](int id) {
    return bridges.ref()[id];
}