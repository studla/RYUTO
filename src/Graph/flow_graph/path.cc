/* 
 * File:   path.cc
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on May 24, 2016, 9:58 AM
 */

#include "path.h"

path::path(exon_edge* e, unsigned int bi, ListDigraph::Arc sa, ListDigraph::Node ln) 
    : identifier(*e), border_index(bi), starting_arc(sa), last_node(ln), last_arc(sa) {

}

path::path(path* p) : identifier(p->identifier), border_index(p->border_index),
        starting_arc(p->starting_arc), last_node(p->last_node), last_arc(p->last_arc){
    
}

path::~path() {
}

