/* 
 * File:   base_manager.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on December 1, 2015, 3:17 PM
 */

#ifndef BASE_MANAGER_H
#define	BASE_MANAGER_H

#include <lemon/list_graph.h>
#include <lemon/lgf_writer.h>
#include <lemon/bfs.h>
#include <set>
#include <unordered_set>
#include "../../Datatype_Templates/maps.h"
#include "../../Datatype_Templates/topsort.h"
#include "../../Datatype_Templates/move.h"
#include "../overlap_graph/overlap_node.h"
#include "../flow_graph/count_raw_node.h"
#include "../flow_graph/count_raw_edge.h"
#include "../flow_graph/region.h"
#include "../flow_graph/edge_types.h"
#include "../flow_graph/edge_length.h"
#include "../pre_graph/pre_graph.h"
#include "../output/alternative_transcript_collection.h"
#include "../output/unsecurity_id.h"
#include "../flow_graph/arc_bridge.h"
#include "../flow_graph/arc_back_bridge.h"
#include "../flow_graph/lpath.h"
#include "../flow_graph/rpath.h"
#include "../flow_graph/capacity_mean.h"
#include "../../Options/options.h"

using namespace lemon;

class base_manager {
public:
    base_manager( pre_graph* raw, exon_meta* meta, const std::string &chromosome);
    virtual ~base_manager();
    
    // filter bins
    void filter_bins();
    // build the splitgraph keeping evidences from multi-exon-atoms
    bool build_extended_splitgraph();
    // build splice graph reducing everything to two-exon splits!
    bool build_basic_splitgraph();
    
    // denoise the graph
    void denoise_graph(ListDigraph::ArcMap<bool> &guided_saves, ListDigraph::ArcMap<bool> &push_block, ListDigraph::ArcMap<bool> &marked_source, ListDigraph::ArcMap<bool> &marked_drain); 
    void denoise_graph_alt(ListDigraph::ArcMap<bool> &guided_saves, ListDigraph::ArcMap<bool> &push_block); 

    void denoise_graph_guided(std::deque<std::deque<ListDigraph::Arc> > &guides);
    
   // void extract_transcripts_from_flow(std::ostream &gs);
    void extract_transcripts_from_flow();
    
    void print_raw( std::deque<overlap_node> &nodes,  std::deque<contained_node> &contained,  std::deque<contained_node> &filtered_contained); 
    
    virtual void print_graph(std::ostream &os);
    virtual void print_transcripts(std::ostream &os);
    
    void print_graph_debug(std::ostream &os, ListDigraph::NodeMap<count_raw_node> &node_counts, ListDigraph::ArcMap<count_raw_edge> &edge_counts);
    
    void print_graph_debug_copy(std::ostream& os, ListDigraph &wc, ListDigraph::ArcMap<capacity_type> &fc, ListDigraph::ArcMap<exon_edge> &ces, ListDigraph::ArcMap<edge_types::edge_type> &cet, ListDigraph::ArcMap<edge_length> &cel, ListDigraph::Node &cs, ListDigraph::Node &ct);
    
    void print_graph_gs(std::ostream& os,
    ListDigraph &wc,
    ListDigraph::ArcMap<capacity_type> &fc,
    ListDigraph::ArcMap<exon_edge> &ces,
    ListDigraph::ArcMap<edge_types::edge_type> &cet,   
    ListDigraph::Node &cs, ListDigraph::Node &ct);
    
    bool has_single_exons();
    void add_single_exons();
    
protected:
    
    // BUILD GRAPH
    
    void update_active( std::deque<overlap_node*> &active_nodes, const unsigned int exon_index);
    void filter_in_graph( std::deque<overlap_node> &nodes, std::deque<contained_node> &contained_nodes);
    void follow_back_erase(overlap_node* on);
    void reduce_transitive( std::deque<overlap_node> &nodes);
    void reduce_transitive_contains( std::deque<contained_node> &contained);
           
    void reduce_single_nodes( std::deque<overlap_node> &nodes,  std::deque<contained_node> &contained, graph_list<exon_group> &raw_exons,
            gmap<exon_group*, overlap_node*> &circle_frag_to_overlap, gmap<exon_group*, contained_node*> &circle_frag_to_contained);
    void follow_single_reduce( overlap_node& current, overlap_node& join,  std::deque<contained_node> &contained,
            gmap<exon_group*, overlap_node*> &circle_frag_to_overlap, gmap<exon_group*, contained_node*> &circle_frag_to_contained);
    
    void create_final_graph( std::deque<overlap_node> &nodes,  std::deque<contained_node> &contained,  std::deque<contained_node> &filtered_contained,
            gmap<exon_group*, overlap_node*> &circle_frag_to_overlap, gmap<exon_group*, contained_node*> &circle_frag_to_contained);
    void add_unaltering_contained_counts( overlap_node* overlap,
         std::deque<contained_node>::iterator &node_it,
         ListDigraph::ArcMap<exon_edge> &edge_specifier,
        ListDigraph::Node* gnodes,
        ListDigraph::ArcMap<count_raw_edge> &edge_counts,
        ListDigraph::NodeMap<count_raw_node> &node_counts);
    void split_edge(unsigned int i,
        range_helper &ranges, // backlinked node that overlaps
        arc_range*& r_it, // range of overlap
        ListDigraph::ArcMap<exon_edge> &edge_specifier,
        ListDigraph::ArcMap<arc_range*> &edge_range,
        ListDigraph::ArcMap<count_raw_edge> &edge_counts,
        ListDigraph::NodeMap<count_raw_node> &node_counts,
        bool* gnode_set, ListDigraph::Node* gnodes, 
         gmap<exon_edge, ListDigraph::Arc> &junction_to_arc,
        std::deque<arc_range> &range_list);
    void split_edge_without_compacting(unsigned int i,
        range_helper &ranges, // backlinked node that overlaps
        arc_range*& r_it, // range of overlap
         ListDigraph::ArcMap<exon_edge> &edge_specifier,
         ListDigraph::ArcMap<arc_range*> &edge_range,
         ListDigraph::ArcMap<count_raw_edge> &edge_counts,
         ListDigraph::NodeMap<count_raw_node> &node_counts,
         ListDigraph::Node* gnodes,
         std::deque<arc_range> &range_list, bool snap_node);
    void add_contained_start(range_helper &ranges,
        unsigned int i,
        ListDigraph::ArcMap<exon_edge> &edge_specifier,
        ListDigraph::NodeMap<count_raw_node> &node_counts,
        ListDigraph::ArcMap<arc_range*> &edge_range,
        ListDigraph::ArcMap<count_raw_edge> &edge_counts,
        std::deque<arc_range> &range_list, bool add_source);
    void add_contained_end(range_helper &ranges,
        unsigned int i,
        ListDigraph::ArcMap<exon_edge> &edge_specifier,
        ListDigraph::NodeMap<count_raw_node> &node_counts,
        ListDigraph::ArcMap<arc_range*> &edge_range,
        ListDigraph::ArcMap<count_raw_edge> &edge_counts,
        std::deque<arc_range> &range_list, bool add_source);
    void snap_contained(range_helper &ranges,
        unsigned int i,
        ListDigraph::ArcMap<exon_edge> &edge_specifier,
        ListDigraph::NodeMap<count_raw_node> &node_counts,
        ListDigraph::ArcMap<arc_range*> &edge_range,
        ListDigraph::ArcMap<count_raw_edge> &edge_counts,
        std::deque<arc_range> &range_list,
        ListDigraph::Node* gnodes, bool* gnode_set,
        gmap<exon_edge, ListDigraph::Arc> &junction_to_arc);
    unsigned int find_index_global_to_sub(unsigned int start, exon_edge* edge);
    
    void compute_edge_length( ListDigraph::Arc arc);
    void create_final_capacities(
        ListDigraph::NodeMap<count_raw_node> &node_counts,
        ListDigraph::ArcMap<count_raw_edge> &edge_counts,
        ListDigraph::Node* gnodes,
        bool* gnode_set);
    
    void enumerate_starts();
    
    void create_region(unsigned int i, rcount startvalue, std::map< rpos,rcount > &lefts, std::map< rpos,rcount > &rights, region &r);
    void create_region_from_node(unsigned int i, count_raw_node &node_counts, region &r);
    void create_region_from_edge(exon_edge& edge, count_raw_edge &edge_counts, region &r);
    
    // FLOW HANDLER
    
    struct flow_paired_path {
        graph_list<paired_exon_group* > frags;
        graph_list<int> in;
        graph_list<int> out;
    };
    struct resolve_count_set {
        
        resolve_count_set() : count(0) {
        };
        rcount count;
        std::set<int> left, right;
    };
    
    struct component {
        component() : edges(0), in_nodes(0), out_nodes(0) {   
        };
        unsigned int edges;
        unsigned int in_nodes;
        unsigned int out_nodes ;
        std::set<int> nodes;
        std::deque< resolve_count_set > fragments;
    };
    
    struct classify_component {
        classify_component() : unassigned_in(0), unassigned_out(0), assigned_in(0), assigned_out(0), assigned_edges(0),
        fully_connected_in(0), fully_connected_out(0), fully_connected_edges(0){  
        };
        unsigned int unassigned_in;
        unsigned int unassigned_out;
        unsigned int assigned_in;
        unsigned int assigned_out;
        unsigned int assigned_edges;
        unsigned int edge_group_offset;
        unsigned int node_group_offset;
        unsigned int fully_connected_in;
        unsigned int fully_connected_out;
        unsigned int fully_connected_edges;
        std::deque<component> components;
    };
    
    struct evidence_group {
        evidence_group() : id(-1), block(false) {};
        evidence_group(int i, bool b) : id(i), block(b) {};
        
        int id;
        bool block; // always set if multiple edges in this group
    };
    
    void find_guide_sources_and_drains(ListDigraph::ArcMap<bool> &guided_saves, std::deque<std::deque<ListDigraph::Arc> > &paths); 
    
    void find_s_t_exons(ListDigraph::ArcMap<bool> &guided_saves, ListDigraph::ArcMap<bool> &marked_source, ListDigraph::ArcMap<bool> &marked_drain, ListDigraph::ArcMap<bool> &neighbouring);
    void filter_introns(ListDigraph::ArcMap<bool> &guided_saves, ListDigraph::ArcMap<bool> &marked_source, ListDigraph::ArcMap<bool> &marked_drain);
    void filter_broken_introns(ListDigraph::ArcMap<bool> &guided_saves, ListDigraph::ArcMap<bool> &marked_source, ListDigraph::ArcMap<bool> &marked_drain);
    
    void prune_sources_and_drains3(ListDigraph::ArcMap<bool> &guided_saves, ListDigraph::ArcMap<bool> &marked_source, ListDigraph::ArcMap<bool> &marked_drain, ListDigraph::ArcMap<bool> &push_block);
    void prune_sources_and_drains4(ListDigraph::ArcMap<bool> &guided_saves, ListDigraph::ArcMap<bool> &marked_source, ListDigraph::ArcMap<bool> &marked_drain, ListDigraph::ArcMap<bool> &push_block);
    void prune_sources_and_drains5(ListDigraph::ArcMap<bool> &guided_saves, ListDigraph::ArcMap<bool> &marked_source, ListDigraph::ArcMap<bool> &marked_drain, ListDigraph::ArcMap<bool> &push_block);

    
    void prune_sources_and_drains_adjacent_starts(ListDigraph::ArcMap<bool> &guided_saves);

    void remove_too_short_source_drain(ListDigraph::ArcMap<bool> &guided_saves, ListDigraph::ArcMap<bool> &marked_source, ListDigraph::ArcMap<bool> &marked_drain, ListDigraph::ArcMap<bool> &consecutive_s_t);
    bool remove_small_low_st(ListDigraph::ArcMap<bool> &guided_saves, ListDigraph::ArcMap<bool> &marked_source, ListDigraph::ArcMap<bool> &marked_drain, ListDigraph::ArcMap<bool> &consecutive_s_t);
    bool erase_low_deviation_st(ListDigraph::ArcMap<bool> &guided_saves, ListDigraph::ArcMap<bool> &marked_source, ListDigraph::ArcMap<bool> &marked_drain, ListDigraph::ArcMap<bool> &consecutive_s_t);
    bool erase_low_boundary_st(ListDigraph::ArcMap<bool> &guided_saves, ListDigraph::ArcMap<bool> &marked_source, ListDigraph::ArcMap<bool> &marked_drain, ListDigraph::ArcMap<bool> &consecutive_s_t);
    bool prune_small_junctions(ListDigraph::ArcMap<bool> &guided_saves, ListDigraph::ArcMap<bool> &marked_source, ListDigraph::ArcMap<bool> &marked_drain, ListDigraph::ArcMap<bool> &consecutive_s_t);
    bool threshold_filter(ListDigraph::ArcMap<bool> &guided_saves, ListDigraph::ArcMap<bool> &marked_source, ListDigraph::ArcMap<bool> &marked_drain, ListDigraph::ArcMap<bool> &consecutive_s_t, float low_mark);

    capacity_type DKMeans( std::vector<capacity_type> &in);
    
    rcount get_forward_median(std::deque<ListDigraph::Node>::iterator &start, std::deque<ListDigraph::Node>::iterator &end,
        ListDigraph::ArcMap<bool> &guided_saves, ListDigraph::ArcMap<bool> &push_block, ListDigraph::ArcMap<bool> &marked_drain);
    rcount get_backward_median(std::deque<ListDigraph::Node>::reverse_iterator &start, std::deque<ListDigraph::Node>::reverse_iterator &end,
        ListDigraph::ArcMap<bool> &guided_saves, ListDigraph::ArcMap<bool> &push_block, ListDigraph::ArcMap<bool> &marked_source);
    
    bool test_consistent_sloping_source(ListDigraph::Node p, ListDigraph::ArcMap<rcount> &maximum, ListDigraph::Arc last, rpos length);
    bool test_consistent_sloping_drain(ListDigraph::Node p, ListDigraph::ArcMap<rcount> &maximum, ListDigraph::Arc last, rpos length);
    
    bool test_consistent_sloping_source(ListDigraph::Node p, rpos length);
    bool test_consistent_sloping_drain(ListDigraph::Node p, rpos length);
    
    bool test_no_pre_sloping_source(ListDigraph::Node p, rpos length, ListDigraph::ArcMap<rcount> &cp_bwd);
    bool test_no_pre_sloping_drain(ListDigraph::Node p, rpos length, ListDigraph::ArcMap<rcount> &cp_bwd);
    
    void prune_unevidenced_arcs(ListDigraph::ArcMap<bool> &guided_saves, ListDigraph::ArcMap<bool> &marked_source, ListDigraph::ArcMap<bool> &marked_drain);
    
    void remove_dead_ends();
    
    void correct_start_end();
    bool recursive_arc_backtracing(exon_edge& goal, unsigned int next_start, unsigned int goal_end_index, ListDigraph::Node n, ListDigraph::Arc la, std::deque<ListDigraph::Arc> &path, bool exon);
    void recursive_start_correction(std::deque< std::pair<std::deque<ListDigraph::Arc>, float> > &start_region, ListDigraph::ArcMap<float > &capacity_correction);
    bool find_starts(ListDigraph::Node p, std::deque< std::pair<std::deque<ListDigraph::Arc>, float> > &start_region, rpos average_fragment_length, std::deque< std::pair<std::deque<ListDigraph::Arc>, float> >::iterator sri);
    void recursive_end_correction(std::deque< std::pair<std::deque<ListDigraph::Arc>, float> > &end_region_region, ListDigraph::ArcMap<float > &capacity_correction);
    bool find_ends(ListDigraph::Node p, std::deque< std::pair<std::deque<ListDigraph::Arc>, float> > &start_region, rpos average_fragment_length, std::deque< std::pair<std::deque<ListDigraph::Arc>, float> >::iterator sri);
    
    void push_potential_forward(ListDigraph::ArcMap<rcount> &cp_fwd, std::deque<ListDigraph::Node> &top_order, ListDigraph::ArcMap<rcount> &average_push, ListDigraph::ArcMap<bool> &push_block);
    void push_potential_backward(ListDigraph::ArcMap<rcount> &cp_fwd, std::deque<ListDigraph::Node> &top_order, ListDigraph::ArcMap<rcount> &average_push, ListDigraph::ArcMap<bool> &push_block);
    void push_potential_forward_alt(ListDigraph::ArcMap<rcount> &cp_fwd, std::deque<ListDigraph::Node> &top_order, ListDigraph::ArcMap<rcount> &average_push, ListDigraph::ArcMap<bool> &push_block);
    void push_potential_backward_alt(ListDigraph::ArcMap<rcount> &cp_fwd, std::deque<ListDigraph::Node> &top_order, ListDigraph::ArcMap<rcount> &average_push, ListDigraph::ArcMap<bool> &push_block);
    float rec_score_forward(ListDigraph::Arc a);
    float rec_score_backward(ListDigraph::Arc a);
    void count_drains_in_distance(ListDigraph::Node p, rpos average_fragment_length, unsigned int &count, unsigned int &calls);
    void count_sources_in_distance(ListDigraph::Node p, rpos average_fragment_length, unsigned int &count, unsigned int &calls);
    void find_over_capacitated_start(ListDigraph::Node p, ListDigraph::ArcMap<rcount> &cp_fwd, std::set<ListDigraph::Node> &nodes);
    void find_over_capacitated_end(ListDigraph::Node p, ListDigraph::ArcMap<rcount> &cp_fwd, std::set<ListDigraph::Node> &nodes);

    
    bool contract_composites(InDegMap<ListDigraph> &in_deg, OutDegMap<ListDigraph> &out_deg);
    void follow_contraction(ListDigraph &wc,
        ListDigraph::ArcMap<capacity_type> &fc, ListDigraph::ArcMap<exon_edge> &ces,
        ListDigraph::ArcMap<edge_types::edge_type> &cet, ListDigraph::ArcMap<edge_length> &cel,
        ListDigraph::ArcMap<capacity_mean> &mc,
        InDegMap<ListDigraph> &in_deg, OutDegMap<ListDigraph> &out_deg,
        ListDigraph::ArcMap<unsigned int> &cycle_id_in, ListDigraph::ArcMap<unsigned int> &cycle_id_out,
        capacity_type path_flow, exon_edge path_ident, edge_length path_length, capacity_mean path_mean,
        unsigned int cidin,
        ListDigraph::Node first, ListDigraph::Node current, bool del) ;
    
    bool simplify_ambiguous(
        ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc,
        ListDigraph::NodeMap<unsecurity_id> &unsecurityId,
        ListDigraph::NodeMap<bool> &resolved,
        ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths,
        InDegMap<ListDigraph> &in_deg, OutDegMap<ListDigraph> &out_deg);
    
    bool simplify_ambiguous_ILP(
        ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc,
        ListDigraph::NodeMap<unsecurity_id> &unsecurityId,
        ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths,
        InDegMap<ListDigraph> &in_deg, OutDegMap<ListDigraph> &out_deg,
        ListDigraph::ArcMap<bool> &barred, bool final_resolve);

  void  compute_ratio(ListDigraph::Node node,
            ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths,
            std::set<int> &component,
            ListDigraph &wc,
            ListDigraph::ArcMap<capacity_type> &fc, float &error, capacity_type &flow);
    
    void find_component_left(std::set<int> &group, int a, ListDigraph &wc, ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths);
    void find_component_right(std::set<int> &group, int a, ListDigraph &wc, ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths);
    
    void extend_path_left(std::deque<path>& paths, path* p, unsigned int& border,
        ListDigraph& wc, ListDigraph::ArcMap<exon_edge>& ces, ListDigraph::ArcMap<edge_types::edge_type> &cet,
	ListDigraph::NodeMap<unsigned int>& cni, ListDigraph::ArcMap<arc_bridge> &know_paths);
    void extend_path_right(std::deque<path>& paths, path* p, unsigned int& border,
        ListDigraph& wc, ListDigraph::ArcMap<exon_edge>& ces, ListDigraph::ArcMap<edge_types::edge_type> &cet,
	ListDigraph::NodeMap<unsigned int>& cni, ListDigraph::ArcMap<arc_bridge> &know_paths);
    
    void compute_edge_groups(std::set<std::set<int> > &hits, std::deque<std::set<int> > &groups);
    
    void resolve_overarching(path_evidence_map<int, path_evidence > &evidences, 
            path_evidence_set<int> &pre_path_match, path_evidence_set<int> &post_path_match,
            path_evidence_map<int, rpos> &post_path_min_length, path_evidence_map<int, rcount> &post_path_count, rpos min_right_hit,
            ListDigraph &wc, ListDigraph::ArcMap<edge_types::edge_type> &cet, ListDigraph::Node n, bool allow_block);
    void resolve_pair(path_evidence_map<int, path_evidence > &evidences, path_evidence_set<int> &pre_path_match,  path_evidence_map<int, rpos> &hit_distance, 
            path_evidence_set<int> &post_path_match, path_evidence_map<int, rpos> &post_path_min_length,
            path_evidence_map<int, std::pair< std::pair< rpos, rpos > ,rpos> > &post_path_range,
            path_evidence_map<int, rcount> &post_path_count,
            ListDigraph &wc, ListDigraph::ArcMap<edge_types::edge_type> &cet, ListDigraph::Node n, bool allow_block);
   
      
public:
    
    static ListDigraph::Arc extract_path(ListDigraph::Arc left_arc, ListDigraph::Arc right_arc, capacity_type capacity,
            ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths,
            ListDigraph &wc,
            ListDigraph::ArcMap<capacity_type> &fc, ListDigraph::ArcMap<capacity_mean> &mc, ListDigraph::ArcMap<exon_edge> &ces,
            ListDigraph::ArcMap<edge_types::edge_type> &cet, ListDigraph::ArcMap<edge_length> &cel,
            ListDigraph::ArcMap<unsigned int> &cycle_id_in, ListDigraph::ArcMap<unsigned int> &cycle_id_out,
            ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc,
            ListDigraph::NodeMap< unsecurity_id > &unsecurityId);
    
    static void unravel_evidences(ListDigraph::Node node,
            ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths,
            ListDigraph &wc,
            ListDigraph::ArcMap<capacity_type> &fc, ListDigraph::ArcMap<capacity_mean> &mc, ListDigraph::ArcMap<exon_edge> &ces,
            ListDigraph::ArcMap<edge_types::edge_type> &cet, ListDigraph::ArcMap<edge_length> &cel,
            ListDigraph::ArcMap<unsigned int> &cycle_id_in, ListDigraph::ArcMap<unsigned int> &cycle_id_out,
            ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc,
            ListDigraph::NodeMap< unsecurity_id > &unsecurityId,
            ListDigraph::NodeMap<bool>* resolved);
    
    static void unravel_evidences_groups(ListDigraph::Node node,
            ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths,
            std::map<int, evidence_group> & left_groups, std::map<int, evidence_group> & right_groups, ListDigraph::ArcMap<bool> &barred,
            ListDigraph &wc,
            ListDigraph::ArcMap<capacity_type> &fc, ListDigraph::ArcMap<capacity_mean> &mc, ListDigraph::ArcMap<exon_edge> &ces,
            ListDigraph::ArcMap<edge_types::edge_type> &cet, ListDigraph::ArcMap<edge_length> &cel,
            ListDigraph::ArcMap<unsigned int> &cycle_id_in, ListDigraph::ArcMap<unsigned int> &cycle_id_out,
            ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc,
            ListDigraph::NodeMap< unsecurity_id > &unsecurityId,  ListDigraph::NodeMap<unsigned int> &ni, unsigned int size);
    
    static capacity_type unravel_evidences_ILP(ListDigraph::Node node,
            ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths,
            std::map<int, evidence_group> & left_groups, std::map<int, evidence_group> & right_groups, ListDigraph::ArcMap<bool> &barred,
            std::deque< resolve_count_set > &hit_counter_all,
            std::set<int> &component,
            ListDigraph &wc,
            ListDigraph::ArcMap<capacity_type> &fc, ListDigraph::ArcMap<capacity_mean> &mc, ListDigraph::ArcMap<exon_edge> &ces,
            ListDigraph::ArcMap<edge_types::edge_type> &cet, ListDigraph::ArcMap<edge_length> &cel,
            ListDigraph::ArcMap<unsigned int> &cycle_id_in, ListDigraph::ArcMap<unsigned int> &cycle_id_out,
            ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc,
            ListDigraph::NodeMap< unsecurity_id> &unsecurityId, ListDigraph::NodeMap<unsigned int> &ni, unsigned int size);
    static void clean_ILP_leftovers(ListDigraph::Node node,
            ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths,
            ListDigraph::ArcMap<bool> &barred, capacity_type barr_limit,
            ListDigraph &wc,
            ListDigraph::ArcMap<capacity_type> &fc, ListDigraph::ArcMap<capacity_mean> &mc, ListDigraph::ArcMap<exon_edge> &ces,
            ListDigraph::ArcMap<edge_types::edge_type> &cet, ListDigraph::ArcMap<edge_length> &cel,
            ListDigraph::ArcMap<unsigned int> &cycle_id_in, ListDigraph::ArcMap<unsigned int> &cycle_id_out,
            ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc,
            ListDigraph::NodeMap< unsecurity_id> &unsecurityId);
    static bool clean_barred_leftovers(ListDigraph::Node node,
            ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths,
            ListDigraph::ArcMap<bool> &barred,
            ListDigraph &wc,
            ListDigraph::ArcMap<capacity_type> &fc, ListDigraph::ArcMap<capacity_mean> &mc, ListDigraph::ArcMap<exon_edge> &ces,
            ListDigraph::ArcMap<edge_types::edge_type> &cet, ListDigraph::ArcMap<edge_length> &cel,
            ListDigraph::ArcMap<unsigned int> &cycle_id_in, ListDigraph::ArcMap<unsigned int> &cycle_id_out,
            ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc,
            ListDigraph::NodeMap< unsecurity_id> &unsecurityId);
    static void unravel_ILP(int left, int right, capacity_type cap_evidence,
            ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths, 
            ListDigraph &wc,
            ListDigraph::ArcMap<capacity_type> &fc, ListDigraph::ArcMap<capacity_mean> &mc, ListDigraph::ArcMap<exon_edge> &ces,
            ListDigraph::ArcMap<edge_types::edge_type> &cet, ListDigraph::ArcMap<edge_length> &cel,
            ListDigraph::ArcMap<unsigned int> &cycle_id_in, ListDigraph::ArcMap<unsigned int> &cycle_id_out,
            ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc,
            ListDigraph::NodeMap< unsecurity_id> &unsecurityId, transcript_unsecurity::evidence evidence);
    static void unravel_single(ListDigraph::Node node,
            ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths,
            ListDigraph &wc,
            ListDigraph::ArcMap<capacity_type> &fc, ListDigraph::ArcMap<capacity_mean> &mc, ListDigraph::ArcMap<exon_edge> &ces,
            ListDigraph::ArcMap<edge_types::edge_type> &cet, ListDigraph::ArcMap<edge_length> &cel,
            ListDigraph::ArcMap<unsigned int> &cycle_id_in, ListDigraph::ArcMap<unsigned int> &cycle_id_out,
            ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc,
            ListDigraph::NodeMap< unsecurity_id> &unsecurityId);
    static void insert_local_known(ListDigraph::Node node, ListDigraph &wc, ListDigraph::ArcMap<arc_bridge> &temp_know_paths, ListDigraph::ArcMap<arc_back_bridge> &temp_know_back_paths,
            ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths);
    static void add_local_known(ListDigraph::Node node, ListDigraph &wc, ListDigraph::ArcMap<arc_bridge> &temp_know_paths, ListDigraph::ArcMap<arc_back_bridge> &temp_know_back_paths,
            ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths);
    static bool bar_smallest_edge(ListDigraph::Node node, ListDigraph &wc, ListDigraph::ArcMap<capacity_type> &fc, ListDigraph::ArcMap<capacity_mean> &mc, ListDigraph::ArcMap<bool> &barred, float ratio, ListDigraph::Arc &barc, float &value, std::unordered_set<int> &block_delete);
    static bool bar_negligible_edge(ListDigraph::Node node, ListDigraph &wc, ListDigraph::ArcMap<capacity_type> &fc, ListDigraph::ArcMap<capacity_mean> &mc, ListDigraph::ArcMap<bool> &barred, float ratio, ListDigraph::Arc &barc, float &value, std::unordered_set<int> &block_delete);
    static bool bar_smallest_edge_by_group(ListDigraph &wc, ListDigraph::ArcMap<capacity_type> &fc, ListDigraph::ArcMap<capacity_mean> &mc, std::set<int> &group, float ratio, ListDigraph::Arc &barc, float &value);

    
    static void unravel_evidence_path_left(int left, int right,
            ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths,
            ListDigraph &wc,
            ListDigraph::ArcMap<capacity_type> &fc, ListDigraph::ArcMap<capacity_mean> &mc, ListDigraph::ArcMap<exon_edge> &ces,
            ListDigraph::ArcMap<edge_types::edge_type> &cet, ListDigraph::ArcMap<edge_length> &cel,
            ListDigraph::ArcMap<unsigned int> &cycle_id_in, ListDigraph::ArcMap<unsigned int> &cycle_id_out,
            ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc,
            ListDigraph::NodeMap< unsecurity_id > &unsecurityId, transcript_unsecurity::evidence evidence,
            ListDigraph::NodeMap<bool>* resolved);
    static void unravel_evidence_path_right(int left, int right,
            ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths, 
            ListDigraph &wc, 
            ListDigraph::ArcMap<capacity_type> &fc, ListDigraph::ArcMap<capacity_mean> &mc, ListDigraph::ArcMap<exon_edge> &ces,
            ListDigraph::ArcMap<edge_types::edge_type> &cet, ListDigraph::ArcMap<edge_length> &cel,
            ListDigraph::ArcMap<unsigned int> &cycle_id_in, ListDigraph::ArcMap<unsigned int> &cycle_id_out,
            ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc,
            ListDigraph::NodeMap< unsecurity_id > &unsecurityId, transcript_unsecurity::evidence evidence,
            ListDigraph::NodeMap<bool>* resolved);
 
    static void unravel_unevidenced_leftovers(ListDigraph::Node &node,
            ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths,
            ListDigraph &wc,
            ListDigraph::ArcMap<capacity_type> &fc, ListDigraph::ArcMap<capacity_mean> &mc, ListDigraph::ArcMap<exon_edge> &ces,
            ListDigraph::ArcMap<edge_types::edge_type> &cet, ListDigraph::ArcMap<edge_length> &cel,
            ListDigraph::ArcMap<unsigned int> &cycle_id_in, ListDigraph::ArcMap<unsigned int> &cycle_id_out,
            ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc,
            ListDigraph::NodeMap< unsecurity_id > &unsecurityId, ListDigraph::ArcMap<bool> *barred);
      
protected:
    
    // VIRTUAL
    
    // virtual function managers MUST implement! 
    virtual void initialize_source_drain_arc(const ListDigraph::Arc &arc) = 0 ;
    virtual const bool expand_exon_nodes() = 0;
    
    virtual void create_node_capacities(ListDigraph::Arc &arc,  region &r) = 0 ;
    virtual void create_edge_capacities(ListDigraph::Arc &arc,  region &r) = 0;
    
    virtual void finalize_flow_graph() = 0;
    
    // creates DAG by breaking cycles and fill flow
    virtual void compute_flow() = 0;
    
    struct single_exon {
        
        single_exon(unsigned int m, capacity_type c, bool g) : meta(m), capacity(c), guide(g) {}
        
        unsigned int meta;
        capacity_type capacity;
        bool guide;
        
        bool operator <(const single_exon& y) const {
            return this->meta < y.meta || this->meta == y.meta && this->capacity < y.capacity;
        }
        
    };
        
    // exon length
    exon_meta* meta;
    
    pre_graph* raw;
    
    std::string chromosome;
    
    // the graph we build up
    ListDigraph g;
    ListDigraph::Node s;
    ListDigraph::Node t;
    
    // we always need those 4 for identifying unresolved nodes
    ListDigraph::ArcMap<exon_edge> edge_specifier; 
    ListDigraph::ArcMap<edge_types::edge_type> edge_type;
    ListDigraph::ArcMap<edge_length> edge_lengths;
    ListDigraph::NodeMap<unsigned int> node_index;
    
    std::set<single_exon> single_exons;
    
    // we have edge capacities for the flow
    ListDigraph::ArcMap<capacity_type> capacity;
    ListDigraph::ArcMap<capacity_mean> means;
    ListDigraph::ArcMap<capacity_type> flow;
    
    // for correction we save up the regions
    ListDigraph::ArcMap<region > regions;
    
    // cycle ids_on edges
    ListDigraph::ArcMap<unsigned int> cycle_id_in;
    ListDigraph::ArcMap<unsigned int> cycle_id_out;
   
public:
    
    // the collection of the output
    alternative_transcript_collection transcripts;
    
};

#endif	/* BASE_MANAGER_H */

