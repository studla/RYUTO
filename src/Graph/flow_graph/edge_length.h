/* 
 * File:   edge_length.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on December 18, 2015, 1:14 PM
 */

#ifndef EDGE_LENGTH_H
#define	EDGE_LENGTH_H

#include <iostream>
#include "../../Datatype_Templates/misc_types.h"

class edge_length {
public:
    edge_length();
    virtual ~edge_length();
    
    rpos first_exon;
    rpos middle;
    rpos last_exon;
    
    rpos without_first();
    rpos without_last();
    
private:

};


std::ostream& operator<<(std::ostream& os, const edge_length& dt);

#endif	/* EDGE_LENGTH_H */

