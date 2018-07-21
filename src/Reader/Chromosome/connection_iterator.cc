/* 
 * File:   connection_iterator.cc
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on February 10, 2016, 4:37 PM
 */

#include "connection_iterator.h"
#include "chromosome.h"
#include "../../Logger/logger.h"

connection_iterator::connection_iterator() : fwd(NULL) , bwd(NULL), in_fwd(true), done(true), counter(0)  {

}


connection_iterator::connection_iterator( chromosome* fwd) : fwd(fwd) , bwd(NULL), in_fwd(true), counter(0) {
   
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("New Iterator " + std::to_string(fwd->chrom_fragments.size())+".\n");
    #endif
    
    if(fwd->chrom_fragments.empty()){
        done = true;
    } else {
        done = false;
        current = fwd->chrom_fragments.begin();
    }
}

connection_iterator::connection_iterator( chromosome* fwd,  chromosome* bwd) : fwd(fwd) , bwd(bwd), counter(0) {
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("New Iterator " + std::to_string(fwd->chrom_fragments.size()) + " " + std::to_string(bwd->chrom_fragments.size())+".\n");
    #endif
    
    if(fwd->chrom_fragments.empty()){
        in_fwd = false;
        if (bwd->chrom_fragments.empty()) {
            done = true;
        } else {
            done = false;
            current = bwd->chrom_fragments.begin();
        }
    } else {
        in_fwd = true;
        done = false;
        current = fwd->chrom_fragments.begin();
        if (bwd->chrom_fragments.empty()) {
            bwd = NULL;
        }
    }
}

connection_iterator::~connection_iterator() {
}


bool connection_iterator::next(greader_list<connected>::iterator &it, unsigned int &index) {
    
    if (done) {
        return false;
    }
    
    index = counter;
    
    if(in_fwd) {
        if (current == fwd->chrom_fragments.end()) {
            if (bwd == NULL) {
                done = false;
                return false;
            } else {
                in_fwd = false;
                current = bwd->chrom_fragments.begin();
            }
        }
    } else {
        if (current == bwd->chrom_fragments.end()) {
            done = true;
            return false;
        }
    }
    
    it = current;
    ++current;
    ++counter;
    
    return true;
}


unsigned int connection_iterator::total() {
    
    unsigned int size = fwd->chrom_fragments.size();
    if (bwd != NULL) {
        size += bwd->chrom_fragments.size();
    }
    
    return size;
}