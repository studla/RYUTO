/* 
 * File:   exon.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on February 3, 2016, 4:13 PM
 */

#ifndef EXON_H
#define	EXON_H

#include "../../Datatype_Templates/misc_types.h"
#include "../../Datatype_Templates/lazy.h"
#include "../../Datatype_Templates/reader_list.h"

class raw_atom;

class exon {
public:
    exon();
    exon(rpos start, rpos end);
    exon(rpos start, rpos end, bool sf, bool ef);
    virtual ~exon();
    
    unsigned int id;
    
    rpos start;
    rpos end;
    
    bool fixed_start;
    bool fixed_end;

    // these are kept for matched reads and do not influence the actual start or end
    //lazy<greader_list<raw_atom*> > matched; // keepng this reference structure would be slower than going without it
    
    bool operator==(const exon &rhs) const {
        return start == rhs.start && end == rhs.end && fixed_start == rhs.fixed_start && fixed_end == rhs.fixed_end;
    }
    
};

bool operator< (const exon& e1, const exon& e2);


#endif	/* EXON_H */

