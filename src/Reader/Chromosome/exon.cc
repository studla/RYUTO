/* 
 * File:   exon.cc
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on February 3, 2016, 4:13 PM
 */

#include <deque>

#include "exon.h"
#include "raw_atom.h"

exon::exon() : id(0), fixed_start(false), fixed_end(false) {
    
}

exon::exon(rpos start, rpos end) : id(0), start(start), end(end), fixed_start(false), fixed_end(false) {
    
}

exon::exon(rpos start, rpos end, bool sf, bool ef) : id(0), start(start), end(end), fixed_start(sf), fixed_end(ef) {
    
}

exon::~exon() {
}


bool operator< (const exon& e1, const exon& e2) {
        
        if (e1.start < e2.start) {
            return true;
        }
        
        if (e1.start > e2.start) {
            return false;
        }
        
        // there should actually be NO overlap, if so, something went wrong
        return e1.end < e2.end;
}