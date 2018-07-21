/* 
 * File:   multi_path_heuristic.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on June 8, 2016, 3:11 PM
 */

#ifndef MULTI_PATH_HEURISTIC_H
#define	MULTI_PATH_HEURISTIC_H

#include "path_finder.h"


class multi_path_heuristic : public path_finder {
public:
    multi_path_heuristic(ListDigraph& wc,
        ListDigraph::Node& s,
        ListDigraph::Node& t,
        ListDigraph::ArcMap<capacity_type>& cfc,
        ListDigraph::ArcMap<capacity_mean> &mc,
        ListDigraph::ArcMap<exon_edge>& ces,
        ListDigraph::ArcMap<edge_types::edge_type>& cet,
        ListDigraph::ArcMap<edge_length> &cel,
        ListDigraph::NodeMap<unsigned int>& cni,
        ListDigraph::ArcMap<arc_bridge>& kp,
        ListDigraph::ArcMap<arc_back_bridge>& kbp,
        ListDigraph::ArcMap<unsigned int>& cycle_id_in,
        ListDigraph::ArcMap<unsigned int>& cycle_id_out, 
        ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc,
        ListDigraph::NodeMap< unsecurity_id> &unsecurityId,
        unsigned int size);
    virtual ~multi_path_heuristic();

    virtual void extract_transcripts( alternative_transcript_collection &results);
    
protected:
    void multi_path_recursion(ListDigraph &wc, ListDigraph::Node &t,
        ListDigraph::ArcMap<capacity_type> &fc, ListDigraph::ArcMap<exon_edge> &ces,
        ListDigraph::ArcMap<edge_types::edge_type> &cet, ListDigraph::ArcMap<edge_length> &cel,
        ListDigraph::NodeMap<unsigned int> &cni,
        ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths,
        ListDigraph::ArcMap<unsigned int>& cycle_id_in, ListDigraph::ArcMap<unsigned int>& cycle_id_out,
        ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc,
        ListDigraph::NodeMap< unsecurity_id > &unsecurityId,
        alternative_transcript_collection &collection,
        alternative_transcript_collection& best, unsigned int* min_count);
    
    void find_min_max( std::deque<std::deque<ListDigraph::Arc> > &path,
        ListDigraph& wc, ListDigraph::Node& s, ListDigraph::ArcMap<capacity_type>& cfc,
        ListDigraph::ArcMap<arc_bridge> &kp);
    

    
    struct evidence_element {
        rcount brigde = 0; // evidence by known paths
        capacity_type cap = 0; // min_max capacity
        std::deque<ListDigraph::Arc> arcs;
        
    };
    
    typedef std::deque<evidence_element> pareto;
    typedef path_evidence_map<int, pareto> dyn_paths;
    
    
    void add_pareto(pareto &pareto_list, evidence_element& insert);
    
};

 

#endif	/* MULTI_PATH_HEURISTIC_H */

