/* 
 * File:   mincost_flow_square.cc
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on January 9, 2017, 11:27 AM
 */

#include "mincost_flow_square.h"

mincost_flow_square::mincost_flow_square(pre_graph* raw, exon_meta* meta, const std::string &chromosome) : mincost_flow_base(raw, meta, chromosome) { 
}


mincost_flow_square::~mincost_flow_square() {
}


void mincost_flow_square::add_offset_edge(capacity_type capacity, capacity_type orig_cap, int exon_count,
        ListDigraph::Node &sn, ListDigraph::Node &tn,
        ListDigraph &og, ListDigraph::ArcMap<capacity_type> &upper, ListDigraph::ArcMap<unsigned_capacity_type> &cost,
        std::deque< ListDigraph::Arc> &reference) {
    
    if (orig_cap == 0) return;
    
    unsigned int mod = 1;
    if (exon_count > 2) {
        mod = exon_count - 2;
    }
    
    capacity_type i = step_interval;
    for (; i <= capacity; i+= step_interval) {
        
        ListDigraph::Arc na = og.addArc(sn, tn);
        upper[na] = step_interval;
        
        capacity_type c1 = i;
        capacity_type c2 = i-step_interval;
        
        cost[na] = (c1*c1 - c2*c2 ) / step_interval * int_scaling_factor / orig_cap * mod;
        reference.push_back(na);
        
    }
    
    i = i - step_interval;
    capacity_type leftover = capacity - i;
    if (leftover > 0) {
        // add one last edge
        ListDigraph::Arc na = og.addArc(sn, tn);
        upper[na] = leftover;
        
        capacity_type c1 = capacity;
        capacity_type c2 = i;
        
        cost[na] = (c1*c1 - c2*c2 ) / leftover * int_scaling_factor / orig_cap * mod;
        reference.push_back(na);
    } 
}

void mincost_flow_square::init_variables(capacity_type max_supply) {
    
    step_interval = max_supply / options::Instance()->get_max_arc_split();
    if (step_interval == 0) {
        step_interval = 1;
    }    
    
     int_scaling_factor = pow(ceil(log10(max_supply)),10);
}