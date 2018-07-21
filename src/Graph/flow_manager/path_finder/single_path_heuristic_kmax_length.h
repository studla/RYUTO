/* 
 * File:   multi_path_heuristic.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on June 8, 2016, 3:11 PM
 */

#ifndef SINGLE_PATH_HEURISTIC_KMAX_LENGTH_H
#define	SINGLE_PATH_HEURISTIC_KMAX_LENGTH_H

#include "path_finder.h"


class single_path_heuristic_kmax_length : public path_finder {
public:
    single_path_heuristic_kmax_length(ListDigraph& wc,
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
        ListDigraph::NodeMap< unsecurity_id > &unsecurityId,
        unsigned int size);
    virtual ~single_path_heuristic_kmax_length();

    virtual void extract_transcripts( alternative_transcript_collection &results);
    
protected:
    
     void find_min_max(std::deque<ListDigraph::Arc> &path);
    

    
    struct evidence_element {
        rcount brigde = 0; // evidence by known paths
        capacity_type cap = 0; // min_max capacity
        rpos length = 0;
        std::deque<ListDigraph::Arc> arcs;
        
    };
    
    typedef std::deque<evidence_element> klist;
    typedef path_evidence_map<int, klist> dyn_paths;
    
    static const unsigned int k = 6;
    
    void add_kmax(klist &pareto_list, evidence_element& insert);
    
};

 

#endif	/* SINGLE_PATH_HEURISTIC_KMAX_LENGTH_H */

