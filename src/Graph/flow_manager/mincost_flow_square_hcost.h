/* 
 * File:   mincost_flow_square.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on January 9, 2017, 11:27 AM
 */

#ifndef MINCOST_FLOW_SQUARE_HCOST_H
#define	MINCOST_FLOW_SQUARE_HCOST_H

#include "mincost_flow_base.h"

class mincost_flow_square_hcost : public mincost_flow_base {
public:
   mincost_flow_square_hcost(pre_graph* raw, exon_meta* meta, const std::string &chromosome, std::set<int> &ids);
   virtual ~mincost_flow_square_hcost();
    
protected:
    
    virtual void add_offset_edge(capacity_type capacity, capacity_type orig_cap, int exon_count,
        ListDigraph::Node &sn, ListDigraph::Node &tn,
        ListDigraph &og, ListDigraph::ArcMap<capacity_type> &upper, ListDigraph::ArcMap<signed_long_capacity_type> &cost,
        std::deque< ListDigraph::Arc> &reference);

    virtual void add_offset_helper(capacity_type capacity,
        ListDigraph::Node &sn, ListDigraph::Node &tn,
        ListDigraph &og, ListDigraph::ArcMap<capacity_type> &upper, ListDigraph::ArcMap<signed_long_capacity_type> &cost,
        std::deque< ListDigraph::Arc> &reference);
    
    virtual void init_variables(capacity_type max_supply);
    
    
    // in this version we create quadratic costs by piecewise approximation, without any caps
    // pseudopolynomial time
    capacity_type step_interval;
    capacity_type int_scaling_factor;
    
};

#endif	/* MINCOST_FLOW_SQUARE_HCOST_H */

