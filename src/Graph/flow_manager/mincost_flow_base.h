/* 
 * File:   mincost_flow_base.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on January 4, 2017, 4:35 PM
 */

#ifndef MINCOST_FLOW_BASE_H
#define	MINCOST_FLOW_BASE_H

#include "base_manager.h"
#include <lemon/cost_scaling.h>
#include <lemon/network_simplex.h>

class mincost_flow_base : public base_manager {
public:
    
    mincost_flow_base(pre_graph* raw, exon_meta* meta, const std::string &chromosome, std::set<int> &ids);
    virtual ~mincost_flow_base();
    
        // interface functions for graph creation
    virtual void initialize_source_drain_arc(const ListDigraph::Arc &arc);
    virtual const bool expand_exon_nodes();
    
    virtual void create_node_capacities(ListDigraph::Arc &arc, flow_series &r);
    virtual void create_edge_capacities(ListDigraph::Arc &arc, flow_series &r);
    
    virtual void finalize_flow_graph(int id);
    
     // interface functions for flow
    virtual void compute_flow(int id);
    
    
    // output with capacity
    virtual void print_graph(std::ostream &os);
    
protected:
    
    // interfaces and functions used for flow computations
    
    virtual void create_initial_offset_graph(
            ListDigraph::NodeMap<unsigned_capacity_type> &ex_flow,
            ListDigraph::NodeMap<ListDigraph::Node > &node_ref,
            ListDigraph::ArcMap<std::deque< ListDigraph::Arc> > &arc_ref_forward,
            ListDigraph::ArcMap<std::deque< ListDigraph::Arc> > &arc_ref_backward,
            capacity_type &max_supply,
            ListDigraph &og,
            ListDigraph::ArcMap<capacity_type> &upper,
            ListDigraph::ArcMap<unsigned_capacity_type> &cost,
            ListDigraph::NodeMap<unsigned_capacity_type> &supply,
            int id);
    
    virtual capacity_type find_exogenous_flow(ListDigraph::NodeMap<unsigned_capacity_type> &ex_flow, int id);
    
    virtual void add_offset_edge(capacity_type capacity, capacity_type orig_cap, int exon_count,
        ListDigraph::Node &sn, ListDigraph::Node &tn,
        ListDigraph &og, ListDigraph::ArcMap<capacity_type> &upper, ListDigraph::ArcMap<unsigned_capacity_type> &cost,
        std::deque< ListDigraph::Arc> &reference);
    
    virtual void add_offset_helper(capacity_type capacity,
        ListDigraph::Node &sn, ListDigraph::Node &tn,
        ListDigraph &og, ListDigraph::ArcMap<capacity_type> &upper, ListDigraph::ArcMap<unsigned_capacity_type> &cost,
        std::deque< ListDigraph::Arc> &reference);
    
    virtual void push_mincost_flow(
            ListDigraph &og,
            ListDigraph::ArcMap<capacity_type> &upper,
            ListDigraph::ArcMap<unsigned_capacity_type> &cost,
            ListDigraph::NodeMap<unsigned_capacity_type> &supply,
            ListDigraph::ArcMap<unsigned_capacity_type> &flowmap);
    
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
    
    virtual void transform_flow_to_orig(ListDigraph::NodeMap<ListDigraph::Node > &node_ref,
            ListDigraph::ArcMap<std::deque< ListDigraph::Arc> > &arc_ref_forward,
            ListDigraph::ArcMap<std::deque< ListDigraph::Arc> > &arc_ref_backward,
            ListDigraph &og,
            ListDigraph::ArcMap<unsigned_capacity_type> &flowmap, int id);
    
    virtual void init_variables(capacity_type max_supply);
    
private:

};

#endif	/* MINCOST_FLOW_BASE_H */

