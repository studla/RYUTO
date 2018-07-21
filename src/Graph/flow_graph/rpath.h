/* 
 * File:   path.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on May 24, 2016, 9:58 AM
 */

#ifndef RPATH_H
#define	RPATH_H

#include <lemon/list_graph.h>
#include "path.h"
#include "edge_length.h"

using namespace lemon;

class rpath : public path {
public:
    
    rpath(exon_edge* e, unsigned int bi, ListDigraph::Arc sa, ListDigraph::Node ln, edge_length& lengths);
    virtual ~rpath();
    
    edge_length lengths;
    
private:

};

#endif	/* RPATH_H */

