/* 
 * File:   overlap_node.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on December 1, 2015, 6:14 PM
 */

#ifndef OVERLAP_NODE_H
#define	OVERLAP_NODE_H

#include <boost/cstdint.hpp>
#include "range_helper.h"
#include "../pre_graph/exon_group.h"

class contained_node;

class overlap_node {
public:
    
    // markers are used for transitive path algorithms
    // also a marker building up the final graph
    enum node_marker : uint8_t { 
        VACANT = 0, INPLAY = 1, ELIMINATED = 2
    };
    
    overlap_node( exon_group* ex);
    virtual ~overlap_node();
    
    void push_link(overlap_node* link);
    void push_backlink(overlap_node* link);
    void erase_backlink(overlap_node* link);
    void erase_link(overlap_node* link);
    void replace_link( overlap_node* old,  overlap_node* nn);
    void replace_backlink( overlap_node* old,  overlap_node* nn);
    // element that is represented by this
    exon_group* exons;
    
    graph_list<contained_node*> contains;
    graph_list<overlap_node*> links;
    graph_list<overlap_node*> back_links; // for overlap reduction
     
    graph_list<bool> edge_marker;
    node_marker marker;
    
    // can be used to disable fragments for including in the graph
    // used for joining nodes into supernode
    bool activated;
    
    // this is a helper object for final graph construction
    // it keeps track of which part is in which graph arc
    range_helper ranges; 
    
private:

    
};

bool operator== (const overlap_node& g1,const  overlap_node& g2);

#endif	/* OVERLAP_NODE_H */

