/* 
 * File:   path.cc
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on May 24, 2016, 9:58 AM
 */

#include "rpath.h"

rpath::rpath(exon_edge* e, unsigned int bi, ListDigraph::Arc sa, ListDigraph::Node ln, edge_length& lengths) 
    : path(e,bi,sa,ln), lengths(lengths) {

}

rpath::~rpath() {
}

