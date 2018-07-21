/* 
 * File:   contained_node.cc
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on December 4, 2015, 10:00 AM
 */

#include <deque>
#include <algorithm>

#include "contained_node.h"
#include "overlap_node.h"

contained_node::contained_node( exon_group* exons) : exons(exons) {
}

contained_node::~contained_node() {
}


void contained_node::add_contained(overlap_node* node) {
    contained_in.push_back(node);
}


void contained_node::replace_contained( overlap_node* old, overlap_node* nc) {

    // remove old elements
    contained_in.erase( std::remove(contained_in.begin(), contained_in.end(), old), contained_in.end());

        // test if new is already in here
    if ( std::find(contained_in.begin(), contained_in.end(), nc) != contained_in.end()) {
        return;
    }
    
    contained_in.push_back(nc);
}


void contained_node::remove_contained( overlap_node* old) {

    // remove old elements
    contained_in.erase( std::remove(contained_in.begin(), contained_in.end(), old), contained_in.end());

}