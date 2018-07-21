/* 
 * File:   exon_meta.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on December 2, 2015, 9:18 AM
 */

#ifndef EXON_META_H
#define	EXON_META_H

#include "../../Datatype_Templates/misc_types.h"
#include <string>
#include "../../Logger/logger.h"

class exon_meta {
public:
    exon_meta();
    virtual ~exon_meta();
    
    void set_size(const unsigned int size);
    
    struct exon_meta_info {
    
        public:

            rpos left, right;
            rpos exon_length;

            // id in the flow graph 
            // start at 0 and keep it tight!
            unsigned int id;

    };
    
    std::string chromosome;
    std::string strand;
    rcount region_reads;
    rcount region_frags;
    
    rcount absolute_reads;
    
    rpos avrg_read_length;
    
    unsigned int size;
    exon_meta_info* exons;
    
    unsigned order_index;
    
    // TODO: add other data needed
    // coverage profile? bias? etc. !
    // possibly GC-Content?
    
    // important! don't set coverage directly in chromosome or reader, 
    // as support for exons might be split up between different multi-hits
    
private:

};

#endif	/* EXON_META_H */

