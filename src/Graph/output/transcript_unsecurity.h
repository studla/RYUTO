/* 
 * File:   transcript_unsecurity.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on May 11, 2016, 1:55 PM
 */

#ifndef TRANSCRIPT_UNSECURITY_H
#define	TRANSCRIPT_UNSECURITY_H

#include <boost/cstdint.hpp>
#include <iostream>
#include "../../Datatype_Templates/misc_types.h"

class transcript_unsecurity {
public:
    
    enum evidence : uint8_t { 
        RESOLVED, EVIDENCED, GUESSED, UNEVIDENCED, BARRED
    };
    
    // RESOLVED: we found prove for this
    // EVIDENCED: unresolved edge but still prove
    // GUESSED: no resolving possible for this node
    // UNEVIDENCED: goes against evidences
    // BARRED: likely false during resolve!
    
    transcript_unsecurity(unsigned int id, evidence e);
    virtual ~transcript_unsecurity();
    
    
    
    unsigned int position; // position of the split (node for exon_meta)
    rpos left;
    rpos right;
    unsigned int id;
    evidence evidenced;
    
private:

};

std::ostream& operator<<(std::ostream& os, const transcript_unsecurity& tu);

bool operator< (const transcript_unsecurity& t1, const transcript_unsecurity& t2);

bool operator== (const transcript_unsecurity& t1, const transcript_unsecurity& t2);

#endif	/* TRANSCRIPT_UNSECURITY_H */

