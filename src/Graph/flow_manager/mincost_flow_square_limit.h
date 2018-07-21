/* 
 * File:   mincost_flow_square_limit.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on January 9, 2017, 2:39 PM
 */

#ifndef MINCOST_FLOW_SQUARE_LIMIT_H
#define	MINCOST_FLOW_SQUARE_LIMIT_H

#include "mincost_flow_square.h"

class mincost_flow_square_limit : public mincost_flow_square {
public:

    mincost_flow_square_limit(pre_graph* raw, exon_meta* meta, const std::string &chromosome);
   virtual ~mincost_flow_square_limit();
    
protected:
    
    virtual void add_offset_edge(capacity_type capacity, capacity_type orig_cap, int exon_count,
        ListDigraph::Node &sn, ListDigraph::Node &tn,
        ListDigraph &og, ListDigraph::ArcMap<capacity_type> &upper, ListDigraph::ArcMap<unsigned_capacity_type> &cost,
        std::deque< ListDigraph::Arc> &reference);
    
    virtual void init_variables(capacity_type max_supply);

    virtual bool modify_repeat(
            ListDigraph::NodeMap<ListDigraph::Node > &node_ref,
            ListDigraph::ArcMap<std::deque< ListDigraph::Arc> > &arc_ref_forward,
            ListDigraph::ArcMap<std::deque< ListDigraph::Arc> > &arc_ref_backward,
            capacity_type &max_supply,
            ListDigraph &og,
            ListDigraph::ArcMap<capacity_type> &upper,
            ListDigraph::ArcMap<unsigned_capacity_type> &cost,
            ListDigraph::NodeMap<unsigned_capacity_type> &supply,
            ListDigraph::ArcMap<unsigned_capacity_type> &flowmap); 
    
    capacity_type overflow_cost;
    capacity_type int_scaling_factor;
};

#endif	/* MINCOST_FLOW_SQUARE_LIMIT_H */

