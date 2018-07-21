/* 
 * File:   path.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on May 24, 2016, 9:58 AM
 */

#ifndef PATH_H
#define	PATH_H

#include <lemon/list_graph.h>
#include "exon_edge.h"

using namespace lemon;

class path {
public:
    
    path(exon_edge* e, unsigned int bi, ListDigraph::Arc sa, ListDigraph::Node ln);
    path(path *p);
    
    virtual ~path();
    
    exon_edge identifier;
    unsigned int border_index;
    ListDigraph::Arc starting_arc;
    ListDigraph::Node last_node;
    ListDigraph::Arc last_arc;
    
private:

};

#endif	/* PATH_H */

