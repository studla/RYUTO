/* 
 * File:   raw_atom.cc
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on February 3, 2016, 4:20 PM
 */

#include <deque>
#include <algorithm>

#include "raw_atom.h"
#include "exon.h"
#include "../../Logger/logger.h"
#include "../../Options/options.h"

raw_atom::raw_atom() : length_filtered(false), source_evidence(false),
        drain_evidence(false), reference_atom(false), has_coverage(false)
{
    length_filtered = options::Instance()->is_min_readlength_set();
}



raw_atom::~raw_atom() {
}

void raw_atom::remove_read(const read_collection* r) {
   reads.ref().erase(std::find(reads.ref().begin(), reads.ref().end(), r));
}


std::string raw_atom::to_string() {
   
    std::ostringstream buffer;
    for( greader_refsorted_list<exon*>::iterator e_it = exons->begin(); e_it != exons->end(); ++e_it) {
        buffer << std::to_string((*e_it)->start) << "-" << std::to_string((*e_it)->end) << " ";
    }
    return buffer.str();
}


// we need to sort atoms unfortunately for compacting them
// since we have not yet any fixation of features, go by exons positions alone
// possibly slow?
bool operator< ( const raw_atom& a1, const  raw_atom& a2) {
     
    greader_refsorted_list<exon*>::const_iterator a1_it = a1.exons.const_ref().begin();
    greader_refsorted_list<exon*>::const_iterator a2_it = a2.exons.const_ref().begin();
    
    for ( ; ; ++a1_it, ++a2_it ) {
    
        if (a1_it == a1.exons.const_ref().end()) { // if both empty without decision -> exons false
            return false;
        }
        
        if (a2_it == a2.exons.const_ref().end()) {
            return true;
        }
        
        // try to move through all cases
        if(*a1_it != *a2_it) {

            // unequal exons, we can decide by postion
            // please note that this is only safe if exons are correct and don't overlap
            // or strictness is violated
            return (*a1_it)->start < (*a2_it)->start;
        }
    
    }  
}

bool operator== (const raw_atom& a1, const raw_atom& a2) {
     
    greader_refsorted_list<exon*>::const_iterator a1_it = a1.exons.const_ref().begin();
    greader_refsorted_list<exon*>::const_iterator a2_it = a2.exons.const_ref().begin();
    
    for ( ; ; ++a1_it, ++a2_it ) {
    
        if (a1_it == a1.exons.const_ref().end() && a2_it == a2.exons.const_ref().end()) { 
            return true;
        }
        
        if(*a1_it != *a2_it) {
            return false;
        }
        
    }  
}