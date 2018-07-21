/* 
 * File:   contained_node.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on December 4, 2015, 10:00 AM
 */

#ifndef CONTAINED_NODE_H
#define	CONTAINED_NODE_H


#include "../pre_graph/exon_group.h"
#include "../../Datatype_Templates/graph_list.h"
#include "range_helper.h"

class overlap_node;

class contained_node {
public:
    contained_node( exon_group* exons);
    virtual ~contained_node();
    
    void add_contained(overlap_node* node);
    void replace_contained( overlap_node* old, overlap_node* nc) ;
    void remove_contained( overlap_node* old) ;
    
    exon_group* exons;
    
    graph_list<overlap_node*> contained_in;
    
    // this is a helper object for final graph construction
    // it keeps track of which part is in which graph arc
    range_helper ranges; 
    
private:

};

#endif	/* CONTAINED_NODE_H */

