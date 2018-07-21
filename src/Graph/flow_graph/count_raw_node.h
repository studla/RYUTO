/* 
 * File:   count_raw_node.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on May 15, 2017, 9:47 AM
 */

#ifndef COUNT_RAW_NODE_H
#define	COUNT_RAW_NODE_H

#include "count_raw_edge.h"

class count_raw_node {
public:
    count_raw_node();

    virtual ~count_raw_node();
    
    std::map< rpos,rcount > lefts;
    rcount total_lefts;
    std::map< rpos,rcount > rights;
    rcount total_rights;
    
    void add_node(count_raw_node* node);
    
    void add_node_start(exon_group* exons);
    void add_node_end(exon_group* exons);
    
    void add_node_index(count_raw_edge& edge, unsigned int i);
    rcount add_node_initial_index(exon_group* exons, unsigned int i, rcount split_v);
    
    capacity_type get_max();
    capacity_type get_min();
    
private:

};

#endif	/* COUNT_RAW_NODE_H */

