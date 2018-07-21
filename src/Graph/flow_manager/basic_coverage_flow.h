/* 
 * File:   basic_coverage_flow.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on April 28, 2016, 2:56 PM
 */

#ifndef BASIC_COVERAGE_FLOW_H
#define	BASIC_COVERAGE_FLOW_H

#include "base_manager.h"
#include <lemon/preflow.h>

class basic_coverage_flow : public base_manager {
public:
    basic_coverage_flow(pre_graph* raw, exon_meta* meta, const std::string &chromosome);
    virtual ~basic_coverage_flow();
    
    // interface functions for graph creation
    virtual void initialize_source_drain_arc(const ListDigraph::Arc &arc);
    virtual const bool expand_exon_nodes();
    
    virtual void create_node_capacities(ListDigraph::Arc &arc, region &r);
    virtual void create_edge_capacities(ListDigraph::Arc &arc, region &r);
    
    virtual void finalize_flow_graph();
    
     // interface functions for flow
    virtual void compute_flow();
    
    // output with capacity
    virtual void print_graph(std::ostream &os);
    
protected:

    // additional methods for resolving cycles
    
    enum circle_marker : uint8_t { 
        NOINFO = 0, SAVE = 1, ELIMINATE = 2
    };
    
    void expand_cycles(ListDigraph &wc, ListDigraph::ArcMap<capacity_type> &fc,
        ListDigraph::ArcMap<capacity_type> &cap, ListDigraph::ArcMap<exon_edge> &ces,
        ListDigraph::ArcMap<edge_types::edge_type> &cet, ListDigraph::ArcMap<edge_length> &cel,
        ListDigraph::NodeMap<unsigned int> &cni,
        ListDigraph::ArcMap<unsigned int> &cycle_id_in, ListDigraph::ArcMap<unsigned int> &cycle_id_out,
        ListDigraph::Node s, ListDigraph::Node t,
        ListDigraph::NodeMap<ListDigraph::Node> &node_ref,
        ListDigraph::ArcMap<ListDigraph::Arc> &arc_ref);
    
    void simplify_ambigous_cycles(ListDigraph &wc, ListDigraph::ArcMap<exon_edge> &ces,
        ListDigraph::ArcMap<edge_types::edge_type> &cet, ListDigraph::NodeMap<unsigned int> &cni,
        std::pair<int,int> arc_pair, 
        graph_list<paired_exon_group* > &evidences,
        gmap< std::pair<int,int> , flow_paired_path > & circles);
    
    void create_leftover_capacities(ListDigraph &wc, ListDigraph::ArcMap<capacity_type> &flow, ListDigraph::ArcMap<capacity_type> &cap, ListDigraph::ArcMap<capacity_type> &new_cap);
    void add_flow(ListDigraph &wc, ListDigraph::ArcMap<capacity_type> &orig_flow, ListDigraph::ArcMap<capacity_type> &add_flow);
    
    void marked_capacities_in(ListDigraph &wc, ListDigraph::ArcMap<capacity_type> &alt_cap,
        std::vector<std::pair<int,int> > &flow_cap, gmap< std::pair<int,int> , flow_paired_path > &circle_links);
    void marked_capacities_out(ListDigraph &wc, ListDigraph::ArcMap<capacity_type> &alt_cap,
        std::vector<std::pair<int,int> > &cap_flow, gmap< std::pair<int,int> , flow_paired_path > &circle_links);
    void marked_capacities_both(ListDigraph &wc, ListDigraph::ArcMap<capacity_type> &alt_cap,
        std::vector<std::pair<int,int> > &cap_cap, gmap< std::pair<int,int> , flow_paired_path > &circle_links);
    
    capacity_type push_circle_flow(ListDigraph &wc, ListDigraph::ArcMap<capacity_type> &alt_cap,
        ListDigraph::ArcMap<capacity_type> &new_flow, capacity_type max_flow, ListDigraph::Node &s, ListDigraph::Node &t);
    
    capacity_type push_graph_flow(ListDigraph &g, ListDigraph::Node s, ListDigraph::Node t, ListDigraph::ArcMap<capacity_type> &cap, ListDigraph::ArcMap<capacity_type> &f) ;
    
    void transfer(ListDigraph &g, ListDigraph::ArcMap<edge_types::edge_type> &cet, ListDigraph::Node s, ListDigraph::Node t, ListDigraph::ArcMap<capacity_type> &cap, ListDigraph::ArcMap<capacity_type> &f, ListDigraph::ArcMap<ListDigraph::Arc> &arc_ref) ;
    
    // unique datastructure
    
    std::deque<ListDigraph::Arc> drain_source_arcs;
    std::set<std::pair<unsigned int, int> > slow_add_sources;
    std::set<std::pair<unsigned int, int> > slow_add_drains;
};

#endif	/* BASIC_COVERAGE_FLOW_H */

