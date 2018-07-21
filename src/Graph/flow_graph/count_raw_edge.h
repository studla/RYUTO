/* 
 * File:   node_count_raw.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on December 15, 2015, 9:59 AM
 */

#ifndef COUNT_RAW_EDGE_H
#define	COUNT_RAW_EDGE_H

#include <vector>
#include "../pre_graph/exon_group.h"
#include "../../Datatype_Templates/misc_types.h"

class count_raw_edge {
public:
    count_raw_edge();
    count_raw_edge(unsigned int size);
    virtual ~count_raw_edge();
      
    unsigned int size; // the length of the edge!
    bool initialized;
    
    // these are spread over the whole of edge
    std::vector<lazy<std::map< rpos,rcount > > > starts;  
    std::vector<lazy<std::map< rpos,rcount > > > ends;
    std::vector<rcount> splits;

    void add_edge(count_raw_edge *edge);
    void add_initial_count(exon_group* exons);
    rcount add_sub_counts_eg(exon_group* exons, int from, int to, rcount split_start);
    rcount add_sub_counts_eg_range(exon_group* exons, int from, int to, unsigned int g1, rcount split_v);
    rcount add_sub_counts_start(exon_group* exons, int g1);
    rcount add_sub_counts_end(exon_group* exons, int g);
    
    void add_sub_counts_re(count_raw_edge &a, int from, int to);
    
    capacity_type get_max();
    capacity_type get_min();
    
private:

};

#endif	/* COUNT_RAW_EDGE_H */

