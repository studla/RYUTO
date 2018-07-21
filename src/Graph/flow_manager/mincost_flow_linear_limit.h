/* 
 * File:   mincost_flow_linear_limit.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on January 11, 2017, 10:30 AM
 */

#ifndef MINCOST_FLOW_LINEAR_LIMIT_H
#define	MINCOST_FLOW_LINEAR_LIMIT_H

#include "mincost_flow_base.h"

class mincost_flow_linear_limit : public mincost_flow_base {
public:

    mincost_flow_linear_limit(pre_graph* raw, exon_meta* meta, const std::string &chromosome);
    virtual ~mincost_flow_linear_limit();
    
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
    
    capacity_type int_scaling_factor;
};

#endif	/* MINCOST_FLOW_LINEAR_LIMIT_H */

