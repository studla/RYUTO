/* 
 * File:   lpath.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on May 30, 2016, 10:16 AM
 */

#ifndef LPATH_H
#define	LPATH_H

#include "path.h"
#include "path_evidence.h"

using namespace lemon;

class lpath : public path {
public:
    lpath(exon_edge* e, unsigned int bi, ListDigraph::Arc sa, ListDigraph::Node ln);
    lpath(lpath* p);
    virtual ~lpath();
    
    void add_evidence(int id, rcount count);
    void add_blocked(int id);
    bool is_blocked(int id);
    
    lpath* parent;
    path_evidence evidence;
    bool is_parent;
    
private:

};

#endif	/* LPATH_H */

