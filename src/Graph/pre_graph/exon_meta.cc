/* 
 * File:   exon_meta.cc
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on December 2, 2015, 9:18 AM
 */

#include "exon_meta.h"

exon_meta::exon_meta(): size(0) {
    
}

void exon_meta::set_size(const unsigned int size)  {
    if (this->size) {
        delete [] exons;
    }
    this->size = size;
    exons = new exon_meta_info[size];
}

exon_meta::~exon_meta() {
    delete [] exons;
}

