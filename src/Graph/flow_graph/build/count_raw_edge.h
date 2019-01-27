/* 
 * File:   node_count_raw.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on December 15, 2015, 9:59 AM
 */

#ifndef COUNT_RAW_EDGE_H
#define	COUNT_RAW_EDGE_H

#include <vector>
#include "../../pre_graph/exon_group.h"
#include "../../../Datatype_Templates/misc_types.h"

class count_raw_edge {
public:
    count_raw_edge();
    count_raw_edge(unsigned int size);
    virtual ~count_raw_edge();
      
    unsigned int size; // the length of the edge! 
    bool initialized;
    
    // these are spread over the whole of edge
    struct series_struct {
        bool initialized = false;
        std::vector<lazy<std::map< rpos,rcount > > > starts;  
        std::vector<lazy<std::map< rpos,rcount > > > ends;
        std::vector<rcount> splits;
    };
    gmap<int, series_struct> series;

    void add_edge(count_raw_edge *edge);
    void add_initial_count(exon_group* exons);
    
    void add_sub_counts_eg(exon_group* exons, int from, int to, gmap<int, rcount> &next_value);
    void add_sub_counts_eg_range(exon_group* exons, int from, int to, unsigned int g1, gmap<int, rcount> &next_value);
    void add_sub_counts_start(exon_group* exons, int g1, gmap<int, rcount> &next_value);
    void add_sub_counts_end(exon_group* exons, int g, gmap<int, rcount> &next_value);
    
    void add_sub_counts_re(count_raw_edge &a, int from, int to);
    
    
private:

};

std::ostream& operator<<(std::ostream& os, const count_raw_edge& cre);

#endif	/* COUNT_RAW_EDGE_H */

