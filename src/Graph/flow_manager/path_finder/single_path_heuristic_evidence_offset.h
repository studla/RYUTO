/* 
 * File:   single_path_heuristic.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on June 8, 2016, 3:03 PM
 */

#ifndef SINGLE_PATH_HEURISTIC_EVIDENCE_OFFSET_H
#define	SINGLE_PATH_HEURISTIC_EVIDENCE_OFFSET_H

#include "path_finder.h"

class single_path_heuristic_evidence_offset : public path_finder {
public:
    single_path_heuristic_evidence_offset(ListDigraph& wc,
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

    virtual ~single_path_heuristic_evidence_offset();
    
    virtual void extract_transcripts( alternative_transcript_collection &results);
    
private:

    void find_min_max(std::deque<ListDigraph::Arc> &path);
};

#endif	/* SINGLE_PATH_HEURISTIC_EVIDENCE_OFFSET_H */

