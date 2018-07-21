/* 
 * File:   connection_iterator.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on February 10, 2016, 4:37 PM
 */

#ifndef CONNECTION_ITERATOR_H
#define	CONNECTION_ITERATOR_H

#include "chromosome.h"

class connection_iterator {
public:
    connection_iterator();
    connection_iterator( chromosome* fwd);
    connection_iterator( chromosome* fwd, chromosome* bwd);
    virtual ~connection_iterator();
    
    chromosome* fwd;
    chromosome* bwd;
    
    bool in_fwd;
    bool done;
    greader_list<connected>::iterator current;
    unsigned int counter;
    
    bool next(greader_list<connected>::iterator &it, unsigned int &index);
    unsigned int total();
    
private:

};

#endif	/* CONNECTION_ITERATOR_H */

