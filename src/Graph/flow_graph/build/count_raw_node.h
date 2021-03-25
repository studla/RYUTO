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
    
    struct series_struct {
        std::map< rpos,rcount > lefts;
        rcount total_lefts = 0;
        std::map< rpos,rcount > rights;
        rcount total_rights = 0;
    };
    gmap<int, series_struct> series;
    
    void add_node(count_raw_node* node);
    
    void add_node_start(exon_group* exons);
    void add_node_end(exon_group* exons);
    
    void add_node_index(count_raw_edge& edge, unsigned int i);
    void add_node_initial_index(exon_group* exons, unsigned int i, gmap<int, rcount> &next_value);
    
private:

};

std::ostream& operator<<(std::ostream& os, const count_raw_node& crn);

#endif	/* COUNT_RAW_NODE_H */

